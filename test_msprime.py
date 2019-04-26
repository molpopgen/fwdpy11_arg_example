import msprime
import numpy as np
import sys
from pathlib import Path
import argparse
from libsequence.msprime import make_SimData
from libsequence.summstats import PolySIM
from libsequence.fst import Fst

def parse_args():
	dstring = "Prototype implementation of ARG tracking and regular garbage collection."
	parser = argparse.ArgumentParser(description=dstring,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--pop1', '-1', nargs=2, default=["tenn","7310"], help="demography type (flat/tenn), initial pop") 
	parser.add_argument('--pop2', '-2', nargs=3, default=[100,110,500], help="size of population 2 in individual diploids, generation after burn-in population 2 arises, generation after burn-in population 2 goes extinct")
	parser.add_argument('--burn_in', '-B', type=int, default=73100.0, help="number of burn-in generations") 
	parser.add_argument('--migration', '-m,', nargs=6, default=[0.1,0.1,(100/14474),111,400,0], help="steady migration rate 1 to 2, steady migration rate 2 to 1, split rate, migration start (after burn-in), migration end (after burn-in), split recovery (population 1 immediately recovers size afer split)") 
	parser.add_argument('--ntheta', '-nT', type=float, default=10.0, help="4Nu: effective mutation rate of neutral mutations scaled to population size 1 at generation 0") 
	parser.add_argument('--selection', '-s', type=float, default=-0.025, help="selection coefficient: -1 < s ") 
	parser.add_argument('--theta', '-T', type=float, default=10.0, help="4Nu: effective mutation rate of selected mutations scaled to population size 1 at generation 0") #for testing against neutral models, set to 0 and let msprime set mutations on the resulting tree
	parser.add_argument('--rho', '-R', type=float, default=10.0, help="4Nr: effective recombination rate scaled to population size 1 at generation 0")
	parser.add_argument('--n_sam1_curr', '-ns1', type=int, default=10, help="Sample size (in diploids) of population 1 in current day.")
	parser.add_argument('--n_sam2_curr', '-ns2', type=int, default=0, help="Sample size (in diploids) of population 2 in current day.")
	parser.add_argument('--anc_sam1', '-as1', nargs='*', default = argparse.SUPPRESS, help="List of ancient samples (generation after burn-in, number of samples - in diploids) of population 1.")
	parser.add_argument('--anc_sam2', '-as2', nargs='*', default = argparse.SUPPRESS, help="List of ancient samples (generation after burn-in, number of samples - in diploids) of population 2.")
	parser.add_argument('--seed', '-S', type=int, default=42, help="RNG seed")
	parser.add_argument('--replicates', '-r', type=int, default=100, help="number of simulation replicates")
	parser.add_argument('--generations', '-g', type=int, default=5920, help="number of generations in flat demography")
	parser.add_argument('--gc', '-G', type=int, default=100, help="GC interval")
	group = parser.add_mutually_exclusive_group(required=False)
	group.add_argument('--init_tree', '-iT', dest='init_tree', action='store_true')
	group.add_argument('--no_init_tree', '-niT', dest='init_tree', action='store_false')
	parser.set_defaults(init_tree=True)
	parser.add_argument('--outfilename', '-o', default="msprime_simulation.txt", help="outfile name")
    
	return parser
	
#burn-in must be 0
def run_msprime(args):
	pop_size_1 = int(args.pop1[1])
	final_generation = args.generations
	split_generation = final_generation - int(args.pop2[1])
	pop2_extinct_gen = final_generation - int(args.pop2[2])
	split_rate = float(args.migration[2])
	split_recovery = bool(float(args.migration[5]))
	pop_size_2 = int(args.pop2[0])
	if(split_rate < 1 and not(split_recovery)):
		pop_size_1 -= pop_size_2
	samples = [msprime.Sample(population=0,time=0)]*(2*args.n_sam1_curr)
	samples.extend([msprime.Sample(population=1,time=0)]*(2*args.n_sam2_curr))
	if(hasattr(args, 'anc_sam1')): 
		i = 0
		while (i < len(args.anc_sam1)):
			samples.extend([msprime.Sample(population=0,time=final_generation-int(args.anc_sam1[i]))]*(2*int(args.anc_sam1[i+1])))
			i += 2
	if(hasattr(args, 'anc_sam2')):
		i = 0
		while (i < len(args.anc_sam2)):
			samples.extend([msprime.Sample(population=1,time=final_generation-int(args.anc_sam2[i]))]*(2*int(args.anc_sam2[i+1])))
			i += 2
	
	population_configurations = []		
	demographic_events = []
	
	if(args.pop1[0] == "flat"):
		population_configurations = [msprime.PopulationConfiguration(initial_size=pop_size_1),msprime.PopulationConfiguration(initial_size=pop_size_2)]
		demographic_events.append(msprime.MassMigration(time=split_generation,source=1,destination=0,proportion=1.0))
		if(split_rate < 1):
			if(not(split_recovery)):
				demographic_events.append(msprime.PopulationParametersChange(time=split_generation, initial_size=(pop_size_1 + pop_size_2), growth_rate=0, population_id=0))
			else:
				demographic_events.append(msprime.PopulationParametersChange(time=split_generation, initial_size=(pop_size_1 - pop_size_2), growth_rate=0, population_id=0))
				demographic_events.append(msprime.PopulationParametersChange(time=(split_generation+1), initial_size=pop_size_1, growth_rate=0, population_id=0))
	else:
		population_configurations = [msprime.PopulationConfiguration(initial_size=512000, growth_rate=0.019552733)]
		demographic_events.append(msprime.PopulationParametersChange(time=205, initial_size=9300, growth_rate=0.003074847, population_id=0))
		demographic_events.append(msprime.PopulationParametersChange(time=920, initial_size=1861, growth_rate=0, population_id=0))
		demographic_events.append(msprime.PopulationParametersChange(time=2040, initial_size=14474, growth_rate=0, population_id=0))
		demographic_events.append(msprime.PopulationParametersChange(time=5920, initial_size=pop_size_1, growth_rate=0, population_id=0))
		
	dd = msprime.DemographyDebugger(population_configurations=population_configurations,demographic_events=demographic_events)
	dd.print_history()
	
	all_sims = msprime.simulate(samples=samples,population_configurations=population_configurations,demographic_events=demographic_events,mutation_rate=args.ntheta/float(4*pop_size_1),recombination_rate=args.rho/float(4*pop_size_1),num_replicates=args.replicates,random_seed=args.seed)
	result_list = []
	
	for sim in all_sims:
		ts_col =  sim.dump_tables()
		
		mysamples = []
		population = [0]
		generation = [0]
		count = 0
		for node in ts_col.nodes:
			if(node.flags == 1):
				if(node.population == population[len(population)-1] and node.time == generation[len(generation)-1]):
					count += 1
				else:
					mysamples.append(count)
					count = 1
					population.append(node.population)
					generation.append(node.time)
			else:
				mysamples.append(count)
				break
		
		pi_array = np.zeros(len(samples))
		fst_array = np.zeros((len(mysamples),len(mysamples)))
		cumsum_samples = np.zeros(1, dtype = np.int64)
		cumsum_samples = np.append(cumsum_samples, np.cumsum(mysamples,dtype=np.int64))

		for i in range(len(mysamples)):
			ts_col = sim.dump_tables()
			sample_nodes = list(range(cumsum_samples[i],cumsum_samples[i+1]))
			ts_col.simplify(samples=sample_nodes)
			subtree_neutral = ts_col.tree_sequence()
				
			sdata = make_SimData(subtree_neutral)
			ps = PolySIM(sdata)
			pi_array[i] = ps.thetapi()/args.ntheta
		
		if(len(mysamples) > 2):
			for i in range(len(mysamples)):
				for j in range((i+1),len(mysamples)):		
					ts_col = sim.dump_tables()
					sample_nodes = list(range(cumsum_samples[i],cumsum_samples[i+1]))
					sample_nodes.extend(list(range(cumsum_samples[j],cumsum_samples[j+1])))
					ts_col.simplify(samples=sample_nodes)
					subtree_neutral = ts_col.tree_sequence()
				
					sdata = make_SimData(subtree_neutral)
					subtree_sample = [mysamples[i],mysamples[j]]
					fst = Fst(sdata,subtree_sample)
					fst_array[i][j] = fst.hsm()/(1-fst.hsm())
					fst_array[j][i] = fst_array[i][j]
				
		elif(len(mysamples) == 2):
			sdata = make_SimData(sim)
			fst = Fst(sdata,mysamples)
			fst_array[0][1] = fst.hsm()/(1-fst.hsm())
			fst_array[1][0] = fst_array[0][1]
		
		result_list.append((fst_array,population,generation,pi_array))
	
	return result_list

if __name__ == "__main__":
	parser = parse_args()
	args = parser.parse_args(sys.argv[1:])
	
	result_list = run_msprime(args)
	
	fst_list = []
	pi_list = []
	for result in result_list:
		pi_list.append(result[3])
		fst_list.append(result[0])
		
	fst_array = np.array(fst_list)
	pi_array = np.array(pi_list)

	mean_fst_array = np.mean(fst_array,axis=0)
	median_fst_array = np.median(fst_array,axis=0)
	std_fst_array = np.std(fst_array,axis=0)
	mean_pi_array = np.mean(pi_array,axis=0)
	median_pi_array = np.median(pi_array,axis=0)
	std_pi_array = np.std(pi_array,axis=0)
	population = result_list[0][1]
	generation = result_list[0][2]
	
	f = open(args.outfilename, "w")
	f.write("population in sample\n")
	np.array(population).tofile(f,sep="\t")
	f.write("\n\ngeneration in sample\n")
	np.array(generation).tofile(f,sep="\t")
	f.write("\n\nexpected_linearlized_fst_array\n")
	for ipop,igen in zip(population,generation):
		for jpop, jgen in zip(population,generation):
			if(ipop == jpop):
				exp_lfst = abs(igen - jgen)/(4*int(args.pop1[1]))
				f.write(str(exp_lfst)+"\t")
			else:
				exp_lfst = (2*(args.generations-int(args.pop2[1]))-igen - jgen)/(4*int(args.pop1[1]))
				f.write(str(exp_lfst)+"\t")	
	
	f.write("\n\nmean_pi_array\n")
	mean_pi_array.tofile(f,sep="\t")
	f.write("\n\nmedian_pi_array\n")
	median_pi_array.tofile(f,sep="\t")
	f.write("\n\nstd_pi_array\n")
	std_pi_array.tofile(f,sep="\t")
	
	f.write("\n\nmean_linearlized_fst_array\n")
	mean_fst_array.tofile(f,sep="\t")
	f.write("\n\nmedian_linearlized_fst_array\n")
	median_fst_array.tofile(f,sep="\t")
	f.write("\n\nstd_linearlized_fst_array\n")
	std_fst_array.tofile(f,sep="\t")

	f.write("\n\nlinearlized_fst_array\t\tpi_array\n")
	for fst_vector, pi_vector in zip(fst_list,pi_list):
		fst_vector.tofile(f,sep="\t")
		f.write("\t")
		pi_vector.tofile(f,sep="\t")
		f.write("\n")

	f.close()
