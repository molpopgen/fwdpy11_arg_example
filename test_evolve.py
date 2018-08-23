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
    n=[init_pop]*(burn_in) #E1: evolve ancestral size to mutation/selection/drift equilibrium
    n.extend([14474]*(5920-2040)) #E2
    n.extend([1861]*(2040-920)) #E3
    n.extend(dem.exponential_size_change(1032,9300,920-205)) #E4
    n.extend(dem.exponential_size_change(9300,512000,205)) #E5
    return np.array(n,dtype=np.uint32)

def parse_args():
    dstring = "Prototype implementation of ARG tracking and regular garbage collection."
    parser = argparse.ArgumentParser(description=dstring,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--pop1', '-1', nargs=3, default=["tenn","7310", "1"], help="demography type (flat/tenn), initial pop, burn-in scale") 
    parser.add_argument('--pop2', '-2', nargs=3, default=[100,110,500], help="size of population 2 in individual diploids, generation after burn-in population 2 arises, generation after burn-in population 2 goes extinct")
    parser.add_argument('--burn_in', '-B', type=int, default=73100.0, help="number of burn-in generations") 
    parser.add_argument('--migration', '-m,', nargs=4,
                        default=[0.1,0.1,111,400], help="migration rate 1 to 2, migration rate 2 to 1, migration start, migration end") 
    parser.add_argument('--ntheta', '-nT', type=float, default=10.0, help="4Nu: effective mutation rate of neutral mutations scaled to population size 1 at generation 0") 
    parser.add_argument('--theta', '-T', type=float, default=10.0, help="4Nu: effective mutation rate of selected mutations scaled to population size 1 at generation 0") #for testing against neutral models, set to 0 and let msprime set mutations on the resulting tree
    parser.add_argument('--rho', '-R', type=float, default=10.0, help="4Nr: effective recombination rate scaled to population size 1 at generation 0")
    parser.add_argument('--n_sam1_curr', '-ns1', type=int, default=10,
                        help="Sample size (in diploids) of population 1 in current day.")
    parser.add_argument('--n_sam2_curr', '-ns2', type=int, default=0,
                        help="Sample size (in diploids) of population 2 in current day.")
    parser.add_argument('--anc_sam1', '-as1', nargs='*', default = argparse.SUPPRESS,
                        help="List of ancient samples (generation, number of samples - in diploids) of population 1.")
    parser.add_argument('--anc_sam2', '-as2', nargs='*', default = argparse.SUPPRESS,
                        help="List of ancient samples (generation, number of samples - in diploids) of population 2.")
    parser.add_argument('--seed', '-S', type=int, default=42, help="RNG seed")
    parser.add_argument('--replicates', '-r', type=int, default=100, help="number of simulation replicates")
    parser.add_argument('--gc', '-G', type=int, default=100, help="GC interval")
	group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--init_tree', '-iT', dest='init_tree', action='store_true')
    group.add_argument('--no_init_tree', '-niT', dest='init_tree', action='store_false')
    parser.set_defaults(init_tree=True)
    
    return parser
    
def run_sim(tuple):
	args = tuple[0]
	seeds = tuple[1]
	init_pop_size = int(args.pop1[1])
	burn_in = args.burn_in
	demography = [init_pop]*(burn_in+5920)
	if(args.pop1[0] == "tenn"):	
	 	demography = get_nlist_tenn(init_pop_size,burn_in)
	 	
	evolver = ea.evolve_track_wrapper(args, demography, seeds)
	print(evolver.times)
	num_sites = evolver.sites.num_rows
	print(num_sites)

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

	msp_rng = msprime.RandomGenerator(seed[3])
	neutral_sites = msprime.SiteTable()
	neutral_mutations = msprime.MutationTable()
	mutgen = msprime.MutationGenerator(msp_rng, args.ntheta/float(4*demography[0])) 
	mutgen.generate(evolver.nodes, evolver.edges, neutral_sites, neutral_mutations)
	num_sites2 = neutral_sites.num_rows
	print(num_sites2)

	trees_neutral = msprime.load_tables(nodes=evolver.nodes, edges=evolver.edges, sites=neutral_sites, mutations=neutral_mutations)

	sdata = make_SimData(trees_neutral)
	
	f = Fst(sd,[5,5])
	
	return f.hsm()
	

if __name__ == "__main__":
	parser = parse_args()
	args = parser.parse_args(sys.argv[1:])
	
	init_pop_size = int(args.pop1[1])
	burn_in = args.burn_in
	args.pop2 = [int(args.pop2[0]),(int(args.pop2[1])+burn_in),(int(args.pop2[2])+burn_in)]
	args.migration = [float(args.migration[0]),float(args.migration[1]),(int(args.migration[2])+burn_in),(int(args.migration[3])+burn_in)]
	
	if(int(args.pop1[1]) <= 0):
		raise RuntimeError("--pop1 initial population size must be > 0")
	if(float(args.pop1[2]) <= 0):
		raise RuntimeError("--pop1 burn-in scale must be > 0")
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
	if(args.migration[0] < 0 or args.migration[0] > 1 or args.migration[1] < 0 or args.migration[1] > 1):
		raise RuntimeError("--migration rates must be between [0,1]")
	if(args.migration[2] > args.migration[3]):
		raise RuntimeError("--migration start must be <= end")
	if(args.pop2[0] > 0 and (args.migration[2] <= args.pop2[1] or args.migration[3] > args.pop2[2] or args.migration[2] > args.pop2[2] or args.migration[3] <= args.pop2[1])):
		raise RuntimeError("--migration start/end must be between pop2 (start,end]")
	if((args.migration[0] > 0 or args.migration[1] > 0) and args.pop2[0] == 0):
		raise RuntimeError("pop2 does not exist, cannot have migration")
	
	# Get 4 seeds for each sim w/0 replacement from [0,1e6)
    	np.random.seed(args.seed)
    	seeds = np.random.choice(range(1000000), 4, replace=False)
	
	tuple = (args,seeds)	
	run_sim(tuple)

