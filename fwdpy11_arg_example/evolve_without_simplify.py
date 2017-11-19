import fwdpy11 as fp11
import fwdpy11.fitness as fp11w
import fwdpy11.model_params as fp11mp
import numpy as np
from .argsimplifier import ArgSimplifier
from .wfarg import AncestryTracker
from .wfarg import evolve_singlepop_regions_track_ancestry
from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions


def evolve(N, rho, simlen=20, seed=42):
    """
    Evolve without selection and without simplifying.

    :param (int) N: Number of diploids
    :param (float) rho: 4Nr
    :param (int) simlen: (Default: 10) Evolve for simlen*N generations.
    :param (int) seed: (Default: 42) RNG seed.

    :rtype: fwdpy11_arg_example.argsimplifier.ArgSimplifier

    :return: Unsimplified data
    """
    pop = fp11.SlocusPop(N)

    # Set up parameters with defaults
    recrate = rho / (4.0 * float(N))
    mutrate_n = 0.0
    mutrate_s = 0.0

    pdict = {'rates': (mutrate_n, mutrate_s, recrate),
            'nregions': [],
             'sregions': [],
             'recregions': [fp11.Region(0, 1, 1)],
             'gvalue': fp11w.SlocusMult(1.0),
             'demography': np.array([N] * simlen * N, dtype=np.uint32)
             }
    params = fp11mp.SlocusParams(**pdict)

    # Set the gc_interval to be so long that it never happens
    gc_interval = 2 * simlen * N
    simplifier = ArgSimplifier(gc_interval, None)
    atracker = AncestryTracker(pop.N, False, 2*pop.N)
    mm = makeMutationRegions(params.nregions, params.sregions)
    rm = makeRecombinationRegions(params.recregions)
    rng = fp11.GSLrng(seed)
    tsim = evolve_singlepop_regions_track_ancestry(rng, pop, atracker, simplifier,
                                                   params.demography,
                                                   params.mutrate_s,
                                                   params.recrate, mm, rm,
                                                   params.gvalue, params.pself)
    # Reverse time:
    atracker.prep_for_gc()
    return atracker
