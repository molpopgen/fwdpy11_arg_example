def evolve_track(rng, pop, params):
    import warnings
    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions
    mm = makeMutationRegions(params.nregions, params.sregions)
    rm = makeRecombinationRegions(params.recregions)

    from .wfarg import evolve_singlepop_regions_track_ancestry
    from .ArgWrapper import ArgWrapper
    ancestry = ArgWrapper()
    evolve_singlepop_regions_track_ancestry(rng, pop, ancestry, 
                                            params.demography,
                                            params.mutrate_s,
                                            params.recrate, mm, rm,
                                            params.gvalue, params.pself)
