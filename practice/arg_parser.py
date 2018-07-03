import argparse
import sys

def parse_args():
    dstring = "Prototype implementation of ARG tracking and regular garbage collection."
    parser = argparse.ArgumentParser(description=dstring,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--pop2', '-2', nargs=3,
                        default=[100,100,200], help="size of population 2 in individual diploids, generation population 2 arises, generation population 2 goes extinct")
    parser.add_argument('--migration', '-m,', nargs=4,
                        default=[0.1,0.1,110,120], help="migration rate 1 to 2, migration rate 2 to 1, migration start, migration end")
    parser.add_argument('--theta', '-T', type=float, default=10.0, help="4Nu")
    parser.add_argument('--rho', '-R', type=float, default=10.0, help="4Nr")
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
	if(args.migration[2] <= args.pop2[1] or args.migration[3] > args.pop2[2] or args.migration[2] > args.pop2[2] or args.migration[3] <= args.pop2[1]):
		raise RuntimeError("--migration start/end must be between pop2 (start,end]")
	if((args.migration[0] > 0 or args.migration[1] > 0) and args.pop2[0] == 0):
		raise RuntimeError("pop2 does not exist, cannot have migration")
	print(args)
    