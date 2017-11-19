#ifndef FWDPY11_ARG_EXAMPLE_EVOLVE_GENERATION_HPP__
#define FWDPY11_ARG_EXAMPLE_EVOLVE_GENERATION_HPP__

#include <cstddef>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <queue>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <fwdpy11/types.hpp>
#include <fwdpy11/rng.hpp>
//#include <fwdpy11/rules/wf_rules.hpp>
#include "ancestry_tracker.hpp"

// This is a copy/paste + modification of fwdpy11's
// existing function to evolve a single-deme, single-region
// population object for one generation.  The modification is to include updating
// an ancestry_tracker object
void
evolve_generation(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    const KTfwd::uint_t N_next, const double mu,
    const KTfwd::traits::mmodel_t<decltype(pop.mutations)>& mmodel,
    const KTfwd::traits::recmodel_t<decltype(pop.gametes),
                                    decltype(pop.mutations)>& recmodel,
    const std::function<std::size_t(const fwdpy11::GSLrng_t&,
                                    const fwdpy11::singlepop_t&)>& pick1,
    const std::function<std::size_t(const fwdpy11::GSLrng_t&,
                                    const fwdpy11::singlepop_t&,
                                    const std::size_t)>& pick2,
    //const update_function& update,
    const std::function<void(const fwdpy11::GSLrng_t&,
                             fwdpy11::singlepop_t::diploid_t&,
                             const fwdpy11::singlepop_t&, const std::size_t,
                             const std::size_t)>& update,
    ancestry_tracker& ancestry);


#endif
