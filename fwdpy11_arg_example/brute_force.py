import fwdpy11
import fwdpy11.fitness
import fwdpy11.model_params
import numpy as np
import msprime
import time
import sys


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
    return atracker


def evolve_track_wrapper(popsize=1000, rho=10000.0, mu=1e-2, seed=42,
                         gc_interval=10,
                         dfe=fwdpy11.ConstantS(0, 1, 1, -0.025, 1.0)):
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
             'demography': np.array([popsize] * 10 * popsize, dtype=np.uint32)
             }

    params = fwdpy11.model_params.SlocusParams(**pdict)
    rng = fwdpy11.GSLrng(seed)
    start_sim = time.time()
    atracker = evolve_track(rng, pop, params, gc_interval)
    stop_sim = time.time()

    sim_nodes = np.array(atracker.nodes, copy=False)
    sim_edges = np.array(atracker.edges, copy=False)
    start_fudge = time.time()
    # Get the node times, convert to float,
    # then convert to backwards in time
    coal_time = np.array(sim_nodes['generation'], dtype=np.float)
    max_gen = coal_time.max()
    samples = [i['id'] for i in sim_nodes if i['generation'] == max_gen]
    coal_time -= max_gen
    coal_time *= -1.0
    stop_fudge = time.time()

    start_msprime = time.time()
    n = msprime.NodeTable()
    n.set_columns(flags=[True for i in range(len(coal_time))],
                  # gives type conversion error from uint32 to int32
                  # without this COPY:
                  population=np.array(sim_nodes['deme'], dtype=np.int32),
                  time=coal_time)
    print("nodes added")
    e = msprime.EdgesetTable()
    for se in sim_edges:
        e.add_row(left=se['left'],
                  right=se['right'],
                  parent=se['parent'],
                  children=(se['child'],))

    # children = [(int(i),) for i in sim_edges['child']]
    # e.set_columns(left=sim_edges['left'],
    #               right=sim_edges['right'],
    #               parent=np.array(sim_edges['parent'], dtype=np.int32),
    #               children=children,
    #               children_length=len(children))
    print("edges added")
    print(n)
    msprime.sort_tables(nodes=n, edgesets=e)
    x = msprime.load_tables(nodes=n, edgesets=e)
    x = x.simplify(sample=samples)
    stop_msprime = time.time()
    return {'sim_time': stop_sim - start_sim,
            'msprime_time': stop_msprime - start_msprime,
            'fudge_time': stop_fudge - start_fudge}
