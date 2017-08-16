#include "evolve_generation.hpp"
#include "ancestry_tracker.hpp"

KTfwd::uint_t
ancestry_recombination_details(
    fwdpy11::singlepop_t& pop, ancestry_tracker& ancestry,
    std::queue<std::size_t>& gamete_recycling_bin,
    const KTfwd::uint_t parental_gamete1, const KTfwd::uint_t parental_gamete2,
	const std::vector<double> & breakpoints,
    const std::tuple<ancestry_tracker::integer_type,
                     ancestry_tracker::integer_type>& pid,
    const std::tuple<ancestry_tracker::integer_type,
                     ancestry_tracker::integer_type>& offspring_indexes)
{
    if (breakpoints.empty())
        {
            ancestry.temp.emplace_back(make_edge(
                0., 1., std::get<0>(pid), std::get<0>(offspring_indexes)));
            return parental_gamete1;
        }
    auto breakpoints_per_parental_chrom = split_breakpoints(breakpoints);
    ancestry.add_edges(breakpoints_per_parental_chrom.first, std::get<0>(pid),
                       std::get<0>(offspring_indexes));

    ancestry.add_edges(breakpoints_per_parental_chrom.second, std::get<1>(pid),
                       std::get<0>(offspring_indexes));
    return KTfwd::recombine_gametes(
        breakpoints, pop.gametes, pop.mutations, parental_gamete1,
        parental_gamete2, gamete_recycling_bin, pop.neutral, pop.selected);
}
