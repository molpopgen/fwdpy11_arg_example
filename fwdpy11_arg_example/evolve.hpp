#ifndef EVOLVE_FUNCTIONS_HPP__
#define EVOLVE_FUNCTIONS_HPP__

#include <fwdpy11/types/SlocusPop.hpp>
#include <pybind11/numpy.h>

double
evolve_track_ancestry(
    const fwdpy11::GSLrng_t& rng, fwdpy11::SlocusPop & pop,
    ancestry_tracker & ancestry, pybind11::function sample_simplify, 
    pybind11::array_t<std::uint32_t> popsizes, 
    pybind11::array_t<std::uint32_t> pop2array, pybind11::array_t<float> migarray,
    const double mu_selected, const double recrate);

#endif
