#ifndef EVOLVE_FUNCTIONS_HPP__
#define EVOLVE_FUNCTIONS_HPP__

#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <pybind11/numpy.h>

double
evolve_singlepop_regions_track_ancestry(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    pybind11::function ancestry_processor,
    pybind11::array_t<std::uint32_t> popsizes, const double mu_selected,
    const double recrate, const KTfwd::extensions::discrete_mut_model& mmodel,
    const KTfwd::extensions::discrete_rec_model& rmodel,
    fwdpy11::single_locus_fitness& fitness, const double selfing_rate);

// double
// evolve_singlepop_regions_track_ancestry_async(
//     const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
//     ancestry_tracker& ancestry, pybind11::function ancestry_processor,
//     const KTfwd::uint_t gc_interval, pybind11::array_t<std::uint32_t> popsizes,
//     const double mu_selected, const double recrate,
//     const KTfwd::extensions::discrete_mut_model& mmodel,
//     const KTfwd::extensions::discrete_rec_model& rmodel,
//     fwdpy11::single_locus_fitness& fitness, const double selfing_rate);
// 
// double
// evolve_singlepop_regions_track_ancestry_python_queue(
//     const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
//     ancestry_tracker& ancestry, pybind11::object python_queue,
//     const KTfwd::uint_t gc_interval, pybind11::array_t<std::uint32_t> popsizes,
//     const double mu_selected, const double recrate,
//     const KTfwd::extensions::discrete_mut_model& mmodel,
//     const KTfwd::extensions::discrete_rec_model& rmodel,
//     fwdpy11::single_locus_fitness& fitness, const double selfing_rate, const int wthreads);

#endif
