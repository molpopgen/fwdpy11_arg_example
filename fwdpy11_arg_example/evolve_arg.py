import fwdpy11
import fwdpy11.fitness
import fwdpy11.model_params
import numpy as np
import msprime


def evolve_track(rng, pop, params, gc_interval):
    import warnings
    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions
    mm = makeMutationRegions(params.nregions, params.sregions)
    rm = makeRecombinationRegions(params.recregions)

    from .wfarg import evolve_singlepop_regions_track_ancestry, AncestryTracker
    from .argsimplifier import ArgSimplifier
    simplifier = ArgSimplifier(gc_interval)
    atracker = AncestryTracker(pop.N)
    evolve_singlepop_regions_track_ancestry(rng, pop, atracker, simplifier,
                                            params.demography,
                                            params.mutrate_s,
                                            params.recrate, mm, rm,
                                            params.gvalue, params.pself)
    if len(atracker.nodes) > 0:
        # TODO
        # The + 1 is b/c we have a bit of a book-keeping
        # thing that we need to document...
        simplifier.simplify(pop.generation+1,atracker)
    return (simplifier,atracker)


def evolve_track_wrapper(popsize=1000, rho=10000.0, mu=1e-2, seed=42,
                         gc_interval=10,
                         dfe=fwdpy11.ConstantS(0, 1, 1, -0.025, 1.0)):
    """
    gc_interval does nothing :)
    """
    if isinstance(dfe, fwdpy11.Sregion) is False:
        raise TypeError("dfe must be a fwdpy11.Sregion")

    if dfe.b != 0.0 or dfe.e != 1.0:
        raise ValueError("DFE beg/end must be 0.0/1.0, repsectively")

    pop = fwdpy11.SlocusPop(popsize)
    recrate = float(rho) / (4.0 * float(popsize))

    pdict = {'rates': (0.0, mu, recrate),
             'nregions': [],
             'sregions': [dfe],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'gvalue': fwdpy11.fitness.SlocusMult(2.0),
             'demography': np.array([popsize] * 20 * popsize, dtype=np.uint32)
             }

    params = fwdpy11.model_params.SlocusParams(**pdict)
    rng = fwdpy11.GSLrng(seed)
    return evolve_track(rng, pop, params, gc_interval)
