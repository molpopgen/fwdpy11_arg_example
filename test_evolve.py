import fwdpy11_arg_example.evolve_arg as ea
import msprime
import numpy as np
import sys
from pathlib import Path
import argparse
import fwdpy11.demography as dem
from libsequence.msprime import make_SimData
from libsequence.fst import Fst
import concurrent.futures

def get_nlist_tenn(init_pop, burn_in):
    """
    Generates a numpy array of the changes in N over time
    There are 5 epochs, with t=0 being the present.
    
    E1: Ne= 7,310  from t=start(8N?) to t = - 5920 generation (Ancestral sizes until 5920 generations ago)
    
    E2: Ne =14,474 from t = -5920 to t = -2040 (Ancient growth at 5920 g ago)
    
    E3: Ne =1,861 from t= -2040 to t= -920 (OOA, first bottle neck 2040 g ago)
    
    E4: Ne = 1,032 to Ne = 9,300 during t = -920 to t = -205 ( second bottle neck and onset of 715 g of exponential growth at rate 0.31% per gen )  
    
    E5: Ne = 9,300 to Ne = 512,000 during t = -205 to t = -0 ( 205 g of exponential growth at rate 1.95% per gen )  
    """
    n=[init_pop]*(burn_in+1) #E1: evolve ancestral size to mutation/selection/drift equilibrium
    n.extend([14474]*(5920-2040)) #E2
    n.extend([1861]*(2040-920)) #E3
    n.extend(dem.exponential_size_change(1032,9300,920-205)) #E4
    n.extend(dem.exponential_size_change(9300,512000,205)) #E5
    return np.array(n,dtype=np.uint32)

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
	parser.add_argument('--outfilename', '-o', default="simulation.txt", help="outfile name")
    
	return parser

def run_sim(tuple):
	args = tuple[0]
	seeds = tuple[1]
	init_pop_size = int(args.pop1[1])
	burn_in = args.burn_in
	demography = [init_pop_size]*(burn_in+args.generations+1)
	if(args.pop1[0] == "tenn"):	
	 	demography = get_nlist_tenn(init_pop_size,burn_in)
	elif(args.migration[2] < 1 and args.migration[5] == False): #for the split recovery model, only the parental generation of population 1 has the bottleneck during the split
		for gen in range(args.pop2[1],len(demography)):
			demography[gen] -=  args.pop2[0]
	
	evolver = ea.evolve_track_wrapper(args, demography, seeds)
	print(evolver.times)
	num_sites = evolver.sites.num_rows
	#print(num_sites)

	# Get a sample of size n_sam1_curr, n_sam2_curr
	final_pop1_size = 2*demography[len(demography)-1]
	final_pop2_size = 0
	if(args.pop2[2] > evolver.pop.generation):
	   final_pop2_size =  args.pop2[0]
	curr_samples = np.random.choice(final_pop1_size, args.n_sam1_curr, replace = False).tolist() #np seed reset in evolve_arg.py 
	if(final_pop2_size > 0 and args.n_sam2_curr > 0):
		curr_samples += (np.random.choice(final_pop2_size, args.n_sam2_curr, replace = False)+final_pop1_size).tolist()
	samples = curr_samples+evolver.anc_samples
	msprime.simplify_tables(samples, nodes = evolver.nodes, edges = evolver.edges, sites = evolver.sites, mutations = evolver.mutations)

	msp_rng = msprime.RandomGenerator(int(seeds[3]))
	neutral_sites = msprime.SiteTable()
	neutral_mutations = msprime.MutationTable()
	mutgen = msprime.MutationGenerator(msp_rng, args.ntheta/float(4*demography[0])) 
	mutgen.generate(evolver.nodes, evolver.edges, neutral_sites, neutral_mutations)
	num_sites2 = neutral_sites.num_rows
	print(num_sites2)
	
	samples = []
	population = [0]
	generation = [0]
	count = 0
	for node in evolver.nodes:
		if(node.flags == 1):
			if(node.population == population[len(population)-1] and node.time == generation[len(generation)-1]):
				count += 1
			else:
				samples.append(count)
				count = 1
				population.append(node.population)
				generation.append(node.time)
		else:
			samples.append(count)
			break

	cumsum_samples = np.zeros(1, dtype = np.int64)
	cumsum_samples = np.append(cumsum_samples, np.cumsum(samples,dtype=np.int64))
	trees_neutral = msprime.load_tables(nodes=evolver.nodes, edges=evolver.edges, sites=neutral_sites, mutations=neutral_mutations)
	
	fst_array = np.zeros((len(samples),len(samples)))
		
	if(len(samples) > 2):
		for i in range(len(samples)):
			for j in range((i+1),len(samples)):
				mynodes = msprime.NodeTable()
				myedges = msprime.EdgeTable()
				mymutations = msprime.MutationTable()
				mysites = msprime.SiteTable()
				
				trees_neutral.dump_tables(nodes=mynodes, edges=myedges, sites=mysites, mutations=mymutations)
				sample_nodes = list(range(cumsum_samples[i],cumsum_samples[i+1]))
				sample_nodes.extend(list(range(cumsum_samples[j],cumsum_samples[j+1])))
				msprime.simplify_tables(samples=sample_nodes, nodes=mynodes, edges=myedges, sites=mysites, mutations=mymutations)
				subtree_neutral = msprime.load_tables(nodes=mynodes, edges=myedges, sites=mysites, mutations=mymutations)
				
				sdata = make_SimData(subtree_neutral)
				subtree_sample = [samples[i],samples[j]]
				fst = Fst(sdata,subtree_sample)
				fst_array[i][j] = fst.hsm()/(1-fst.hsm())
				fst_array[j][i] = fst_array[i][j]
				
	elif(len(samples) == 2):
	
		sdata = make_SimData(trees_neutral)
		fst = Fst(sdata,samples)
		fst_array[0][1] = fst.hsm()/(1-fst.hsm())
		fst_array[1][0] = fst_array[0][1]
		
	return (fst_array,population,generation)		

if __name__ == "__main__":
	parser = parse_args()
	args = parser.parse_args(sys.argv[1:])
	
	init_pop_size = int(args.pop1[1])
	burn_in = args.burn_in
	args.pop2 = [int(args.pop2[0]),(int(args.pop2[1])+burn_in),(int(args.pop2[2])+burn_in)]
	split_rate = float(args.migration[2])
	if(split_rate != 1 and split_rate != 0 and args.pop1[0] == "flat"):
		split_rate = min(args.pop2[0]/init_pop_size,1.0)
	if(args.pop1[0] == "tenn" and split_rate != 0):
		split_rate = 1
			
	args.migration = [float(args.migration[0]),float(args.migration[1]),split_rate,(int(args.migration[3])+burn_in),(int(args.migration[4])+burn_in),bool(float(args.migration[5]))]
	
	if(hasattr(args, 'anc_sam1')): 
		args.anc_sam1 = [int(i)+burn_in*(j%2==0) for j,i in enumerate(args.anc_sam1)] #add burn-in generation to sample generations
	if(hasattr(args, 'anc_sam2')):
		args.anc_sam2 = [int(i)+burn_in*(j%2==0) for j,i in enumerate(args.anc_sam2)]
	if(args.selection <= -1):
		raise RuntimeError("--selection coefficient must be > -1")
	if(int(args.pop1[1]) <= 0):
		raise RuntimeError("--pop1 initial population size must be > 0")
	if(args.pop2[0] < 0):
		raise RuntimeError("--pop2 pop_size must be >= 0")
	if(args.pop2[0] > 0 and args.pop2[1] < 0):
		raise RuntimeError("--pop2 must arise after generation 0")
	if(args.pop2[1] > args.pop2[2]):
		raise RuntimeError("--pop2 origin must be <= extinction")
	if(hasattr(args, 'anc_sam1') and len(args.anc_sam1) % 2):
		raise RuntimeError("--anc_sam1 must have the generation and the number of samples taken")
	if(hasattr(args, 'anc_sam2') and len(args.anc_sam2) % 2):
		raise RuntimeError("--anc_sam2 must have the generation and the number of samples taken")
	if(args.migration[0] < 0 or args.migration[0] > 1 or args.migration[1] < 0 or args.migration[1] > 1 or args.migration[2] < 0 or args.migration[2] > 1 or (args.migration[2] == 0 and args.pop2[0] > 0)):
		raise RuntimeError("--steady migration rates must be between [0,1] and initial migration rate to population 2 must be (0,1] if population 2 exists")
	if(args.migration[3] > args.migration[4]):
		raise RuntimeError("--migration start must be <= end")
	if((args.migration[0] > 0 or args.migration[1] > 0 or args.migration[2] > 0) and args.pop2[0] == 0):
		raise RuntimeError("pop2 does not exist, cannot have migration")
	if(args.migration[0] > 0 and args.migration[1] > 0 and (args.migration[3] <= args.pop2[1] or args.migration[4] > args.pop2[2] or args.migration[3] > args.pop2[2] or args.migration[4] <= args.pop2[1])):
		raise RuntimeError("--migration start/end must be between pop2 (start,end]")
	if(args.replicates <= 0):
		raise RuntimeError("number of replicates must be >= 1")
	# Get 4 seeds for each sim w/0 replacement from [0,1e6)
	np.random.seed(args.seed)
	seeds = np.random.choice(range(1000000), 4*args.replicates, replace=False)

	seed_list = [(seeds[i],seeds[i+1],seeds[i+2],seeds[i+3]) for i in range(0,len(seeds),4)]
	#run_sim((args,seed_list[0]))
	
	result_list = []
	with concurrent.futures.ProcessPoolExecutor() as pool:
		futures = {pool.submit(run_sim, (args,i)) for i in seed_list}
		for fut in concurrent.futures.as_completed(futures):
			result_list.append(fut.result())

	fst_list = []
	for result in result_list:
		fst_list.append(result[0])
		
	fst_array = np.array(fst_list)
	mean_fst_array = np.mean(fst_array,axis=0)
	median_fst_array = np.median(fst_array,axis=0)
	std_fst_array = np.std(fst_array,axis=0)
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
				exp_lfst = (2*(args.generations-int(args.pop2[1])+burn_in)-igen - jgen)/(4*int(args.pop1[1]))
				f.write(str(exp_lfst)+"\t")
	f.write("\n\nmean_linearlized_fst_array\n")
	mean_fst_array.tofile(f,sep="\t")
	f.write("\n\nmedian_linearlized_fst_array\n")
	median_fst_array.tofile(f,sep="\t")
	f.write("\n\nstd_linearlized_fst_array\n")
	std_fst_array.tofile(f,sep="\t")
	f.write("\n\nlinearlized_fst_array\n")
	for fst_array in fst_list:
		fst_array.tofile(f,sep="\t")
		f.write("\n")
	f.close()
