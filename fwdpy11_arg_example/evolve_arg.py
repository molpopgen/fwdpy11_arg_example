import fwdpy11
import fwdpy11.model_params
import fwdpy11.genetic_values
import numpy as np
import msprime

class sampler(object):
    def __init__(self, samples_pop1, samples_pop2, seed):
        np.random.seed(seed)
        self.__samples_pop1 = samples_pop1
        self.__samples_pop2 = samples_pop2
        self.__samples_index_pop1 = 0
        self.__samples_index_pop2 = 0

    def __call__(self, generation, pop_size1, pop_size2, params, total_generations):
        samples = np.array([],dtype=np.int64)
        if(self.__samples_index_pop1+1 < len(self.__samples_pop1) and generation == self.__samples_pop1[self.__samples_index_pop1]):
        	samples = np.append(samples,np.random.choice(int(pop_size1), self.__samples_pop1[self.__samples_index_pop1+1], replace=False))
        	self.__samples_index_pop1 += 2
        if(self.__samples_index_pop2+1 < len(self.__samples_pop2) and generation == self.__samples_pop2[self.__samples_index_pop2]):
        	samples = np.append(samples,(np.random.choice(int(pop_size2), self.__samples_pop2[self.__samples_index_pop2+1], replace=False)+int(pop_size1)))
        	self.__samples_index_pop2 += 2
        return samples

def evolve_track(rng, args, pop, params, seeds, init_with_TreeSequence):
    """
    Evolve a population and track its ancestry using msprime.

    :param rng: A fwdpy11.GSLrng
    :param pop: A fwdpy11.SlocusPop
    :param params: A fwdpy11.SlocusParams
    :param gc_interval: An integer representing how often to simplify the ancestry.

    :rtype: tuple

    :return: An instance of ARGsimplifier.
    """
    import warnings
    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    # Enforce min left end of 0.0, which is an msprime
    # requirement:
    if any(i.b < 0.0 for i in params.nregions) is True:
        raise RuntimeError("Minimum possible position is 0.0")

    if any(i.b < 0.0 for i in params.recregions) is True:
        raise RuntimeError("Minimum possible position is 0.0")

    if any(i.b < 0.0 for i in params.sregions) is True:
        raise RuntimeError("Minimum possible position is 0.0")

    from .argevolver import ArgEvolver
    initial_TreeSequence = None
       
    if init_with_TreeSequence is True:
        initial_TreeSequence = msprime.simulate(
            2 * pop.N, recombination_rate=params.recrate / 2.0, Ne=pop.N, random_seed=seeds[1])
    
    samples_pop1 =[] 
    samples_pop2 =[] 
    if(hasattr(args, 'anc_sam1')): 
        samples_pop1 = args.anc_sam1
    if(hasattr(args, 'anc_sam2')):
        samples_pop2 = args.anc_sam2
    return ArgEvolver(rng, args, pop, params, sampler(samples_pop1,samples_pop2,seeds[2]), initial_TreeSequence)


def evolve_track_wrapper(args, demography, seeds):
    """
    Wrapper around evolve_track to facilitate testing.

    :param popsize: Diploid population size.
    :param rho: 4Nr
    :param mu: Mutation rate to selected alleles
    :param seed: RNG seed
    :param gc_interval: Garbage collection interval.
    :param dfe: An instance of a fwdpy11.Sregion

    :rtype: tuple

    :return: See evolve_track for details.
    """
    
    initial_popsize = demography[0]
    pop = fwdpy11.SlocusPop(initial_popsize)
    recrate = float(args.rho) / (4.0 * float(initial_popsize))
    mu = float(args.theta) / (4.0 * float(initial_popsize))
    
    dfe = [fwdpy11.ConstantS(0, 1, 1, args.selection, 1.0)]
    
    if(not(args.single_locus)):
    	dfe = [fwdpy11.ConstantS(0, args.region_breaks[0], 1, args.selection, 1.0),
    		   fwdpy11.ConstantS(args.region_breaks[1], 1, 1, args.selection, 1.0)]
    
    if isinstance(dfe[0], fwdpy11.Sregion) is False:
        raise TypeError("dfe must be a fwdpy11.Sregion")
    if(not(args.single_locus)):
        if isinstance(dfe[1], fwdpy11.Sregion) is False:
        	raise TypeError("dfe must be a fwdpy11.Sregion")
        	
    recregion = [fwdpy11.Region(0, 1, 1)]
    if(not(args.single_locus)):
        recregion = [fwdpy11.Region(0, args.region_breaks[0], 1),
        			 fwdpy11.Region(args.region_breaks[1], 1, 1)]
    
    pdict = {'nregions': [],
                'sregions': dfe,
                'recregions': recregion,
                'rates': (0.0, mu, recrate),
                'demography': demography,
                'gvalue': fwdpy11.genetic_values.SlocusMult(2.0) #this value does not get used, this is set in evolve.cc: const auto ff = fwdpp::multuplicative_diploid(fwdpp::fitness(2.0)); 
            }
            
    params = fwdpy11.model_params.ModelParams(**pdict)
    rng = fwdpy11.GSLrng(seeds[0])
    
    return evolve_track(rng, args, pop, params, seeds, args.init_tree)
