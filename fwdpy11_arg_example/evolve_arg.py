import fwdpy11
import fwdpy11.fitness
import fwdpy11.model_params
import numpy as np
import msprime


def evolve_track(rng, pop, params, gc_interval, init_with_TreeSequence=False, msprime_seed=None, async=False, use_queue=False, qsize=2):
    """
    Evolve a population and track its ancestry using msprime.

    :param rng: A fwdpy11.GSLrng
    :param pop: A fwdpy11.SlocusPop
    :param params: A fwdpy11.SlocusParams
    :param gc_interval: An integer representing how often to simplify the ancestry.

    :rtype: tuple

    :return: An instance of ARGsimplifier, an instance of AncestryTracker, and the total time spent simulating.
    """
    import warnings
    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    if async is True and use_queue is True:
        raise ValueError("async and use_queue cannot both be True")

    # Enforce min left end of 0.0, which is an msprime
    # requirement:
    if any(i.b < 0.0 for i in params.nregions) is True:
        raise RuntimeError("Minimum possible position is 0.0")

    if any(i.b < 0.0 for i in params.recregions) is True:
        raise RuntimeError("Minimum possible position is 0.0")

    if any(i.b < 0.0 for i in params.sregions) is True:
        raise RuntimeError("Minimum possible position is 0.0")

    from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions
    mm = makeMutationRegions(params.nregions, params.sregions)
    rm = makeRecombinationRegions(params.recregions)

    from .wfarg import evolve_singlepop_regions_track_ancestry, AncestryTracker
    from .wfarg import evolve_singlepop_regions_track_ancestry_async
    from .wfarg import evolve_singlepop_regions_track_ancestry_python_queue
    from .argsimplifier import ArgSimplifier
    initial_TreeSequence = None
    next_index = 2 * pop.N
    if init_with_TreeSequence is True:
        if msprime_seed is None:
            import warnings
            warnings.warn(
                "msprime_seed is None. Results will not be reprodicible.")
        initial_TreeSequence = msprime.simulate(
            2 * pop.N, recombination_rate=params.recrate / 2.0, Ne=pop.N, random_seed=msprime_seed)
        next_index = initial_TreeSequence.num_nodes
    simplifier = ArgSimplifier(gc_interval, initial_TreeSequence)
    atracker = AncestryTracker(pop.N, init_with_TreeSequence, next_index)
    if async is False and use_queue is False:
        tsim = evolve_singlepop_regions_track_ancestry(rng, pop, atracker, simplifier,
                                                       params.demography,
                                                       params.mutrate_s,
                                                       params.recrate, mm, rm,
                                                       params.gvalue, params.pself)
    elif async is True:
        tsim = evolve_singlepop_regions_track_ancestry_async(rng, pop, atracker, simplifier,
                                                             gc_interval,
                                                             params.demography,
                                                             params.mutrate_s,
                                                             params.recrate, mm, rm,
                                                             params.gvalue, params.pself)
    else:
        import threading
        import queue
        q = queue.Queue(4)

        def worker():
            while True:
                data = q.get()
                if data is None:
                    break;
                print(type(data), data[0])
                simplifier(data[0], data[1])
                q.task_done()
        t = threading.Thread(target=worker)
        t.start()
        tsim = evolve_singlepop_regions_track_ancestry_python_queue(rng, pop, atracker, q,
                                                             gc_interval,
                                                             params.demography, params.mutrate_s,
                                                             params.recrate, mm, rm, params.gvalue, params.pself)
        q.put(None)
        t.join()

    print("we have returned from the C++ side of the sim")
    if len(atracker.nodes) > 0:
        print("final simplify")
        # TODO
        # The + 1 is b/c we have a bit of a book-keeping
        # thing that we need to document...
        simplifier.simplify(pop.generation + 1, atracker)
    return (simplifier, atracker, tsim)


def evolve_track_wrapper(popsize=1000, rho=10000.0, mu=1e-2, seed=42,
                         gc_interval=10,
                         dfe=fwdpy11.ConstantS(0, 1, 1, -0.025, 1.0)):
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
    if isinstance(dfe, fwdpy11.Sregion) is False:
        raise TypeError("dfe must be a fwdpy11.Sregion")

    if dfe.b != 0.0 or dfe.e != 1.0:
        raise ValueError("DFE beg/end must be 0.0/1.0, repsectively")

    pop=fwdpy11.SlocusPop(popsize)
    recrate=float(rho) / (4.0 * float(popsize))

    pdict={'rates': (0.0, mu, recrate),
             'nregions': [],
             'sregions': [dfe],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'gvalue': fwdpy11.fitness.SlocusMult(2.0),
             'demography': np.array([popsize] * 20 * popsize, dtype=np.uint32)
             }

    params=fwdpy11.model_params.SlocusParams(**pdict)
    rng=fwdpy11.GSLrng(seed)
    return evolve_track(rng, pop, params, gc_interval)
