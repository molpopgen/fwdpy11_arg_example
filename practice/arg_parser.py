import argparse
import sys
import fwdpy11.demography as dem
import numpy as np

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
    parser.add_argument('--pop1', '-1', nargs=3,
                        default=["tenn","7310", "1"], help="demography type (flat/tenn), initial pop, burn-in scale") 
    parser.add_argument('--pop2', '-2', nargs=3,
                        default=[100,110,500], help="size of population 2 in individual diploids, generation after burn-in population 2 arises, generation after burn-in population 2 goes extinct") 
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
    parser.add_argument('--gc', '-G', type=int,
                        default=100, help="GC interval")

    return parser

    
    
if __name__ == "__main__":
	parser = parse_args()
	args = parser.parse_args(sys.argv[1:])
	print(args)
	
	init_pop_size = int(args.pop1[1])
	burn_in = int(float(args.pop1[2])*init_pop_size)
	args.pop2 = [int(args.pop2[0]),(int(args.pop2[1])+burn_in),(int(args.pop2[2])+burn_in)]
	args.migration = [float(args.migration[0]),float(args.migration[1]),(int(args.pop2[1])+burn_in),(int(args.pop2[2])+burn_in)]
	
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
	
	print(args)
	
	
	demography = [init_pop_size]*(burn_in)
	if(args.pop1[0] == "tenn"):	
	 	demography = get_nlist_tenn(init_pop_size,burn_in)
	print(demography)
	print(args)
    