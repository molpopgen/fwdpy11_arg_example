import fwdpy11
import fwdpy11.fitness
import fwdpy11.model_params
import numpy as np
import msprime

class sampler(object):
    def __init__(self, samples_pop1, samples_pop2, seed = None):
        np.random.seed(seed)
        if seed is None:
           import warnings
           warnings.warn("sampler seed is None. Results will not be reprodicible.")
        self.__samples_pop1 = samples_pop1
        self.__samples_pop2 = samples_pop2
        self.__samples_index_pop1 = 0
        self.__samples_index_pop2 = 0

    def __call__(self, generation, pop_size1, pop_size2, params, total_generations):
        samples = np.array([])
        if(self.__samples_index_pop1+1 < len(self.__samples_pop1) and generation == self.__samples_pop1[self.__samples_index_pop1]):
        	self.__samples_index_pop1 += 1
        	np.append(samples,np.random.choice(int(pop_size1), self.__samples_pop1[self.__samples_index_pop1+1], replace=False))
        if(self.__samples_index_pop2+1 < len(self.__samples_pop2) and generation == self.__samples_pop2[self.__samples_index_pop2]):
        	self.__samples_index_pop2 += 1
        	np.append(samples,(np.random.choice(int(pop_size2), self.__samples_pop2[self.__samples_index_pop2+1], replace=False)+int(pop_size1)))
        return samples

def evolve_track(rng, parsed_args, pop, params, seed=None, init_with_TreeSequence=False, msprime_seed=None):
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
        if msprime_seed is None:
            import warnings
            warnings.warn(
                "msprime_seed is None. Results will not be reprodicible.")
        initial_TreeSequence = msprime.simulate(
            2 * pop.N, recombination_rate=params.recrate / 2.0, Ne=pop.N, random_seed=msprime_seed)
    
    samples_pop1 =[]
    samples_pop2 =[]
    if(hasattr(parsed_args, 'anc_sam1')): 
        samples_pop1 = parsed_args.anc_sam1
    if(hasattr(parsed_args, 'anc_sam2')):
        samples_pop2 = parsed_args.anc_sam2
    return ArgEvolver(rng, parsed_args, pop, params, sampler(samples_pop1,samples_pop2,seed), initial_TreeSequence)


def evolve_track_wrapper(parsed_args, demography):
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
    
    dfe = fwdpy11.ConstantS(0, 1, 1, -0.025, 1.0)
    
    if isinstance(dfe, fwdpy11.Sregion) is False:
        raise TypeError("dfe must be a fwdpy11.Sregion")

    if dfe.b != 0.0 or dfe.e != 1.0:
        raise ValueError("DFE beg/end must be 0.0/1.0, repsectively")

    initial_popsize = demography[0]
    pop = fwdpy11.SlocusPop(initial_popsize)
    recrate = float(parsed_args.rho) / (4.0 * float(initial_popsize))
    mu = float(parsed_args.theta) / (4.0 * float(initial_popsize))
    seed = parsed_args.seed
	
    pdict = {'rates': (0.0, mu, recrate),
             'nregions': [],
             'sregions': [dfe],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'gvalue': fwdpy11.fitness.SlocusMult(2.0),
             'demography': demography
             }

    params = fwdpy11.model_params.SlocusParams(**pdict)
    rng = fwdpy11.GSLrng(seed)
    return evolve_track(rng, parsed_args, pop, params, seed, True, seed+3)
