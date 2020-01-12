#ifndef FWDPY11_ARG_EXAMPLE_HANDLE_RECOMBINATION_HPP__
#define FWDPY11_ARG_EXAMPLE_HANDLE_RECOMBINATION_HPP__

#include <tuple>
#include <vector>
#include "ancestry_tracker.hpp"

std::pair<std::vector<std::pair<double, double>>,
          std::vector<std::pair<double, double>>>
split_breakpoints(const std::vector<double>& breakpoints,
                  const double start = 0., const double stop = 1.);

// breakpoints is passed in non-const.  We filter out
// double x-overs in this function. Really bad things will
// happen if breakpoints is not up to fwdpp's spec!
fwdpp::uint_t ancestry_rec_mut_details(
    fwdpy11::SlocusPop& pop, ancestry_tracker& ancestry,
    std::queue<std::size_t>& gamete_recycling_bin,
    const fwdpp::uint_t parental_gamete1, const fwdpp::uint_t parental_gamete2,
    std::vector<double>& breakpoints,
    const std::vector<fwdpp::uint_t>& new_mutations,
    const std::tuple<ancestry_tracker::integer_type,
                     ancestry_tracker::integer_type>& pid,
    const ancestry_tracker::integer_type offspring_index);

#endif
