#include "handle_mut_rec.hpp"
#include <fwdpp/debug.hpp>
#include <fwdpp/mutate_recombine.hpp>

namespace py = pybind11;

std::pair<std::vector<std::pair<double, double>>,
          std::vector<std::pair<double, double>>>
split_breakpoints(const std::vector<double>& breakpoints, const double start,
                  const double stop)
{
    using break_array = std::vector<std::pair<double, double>>;
    if (breakpoints.empty())
    	return std::make_pair(break_array{std::make_pair(0.,1.)},break_array{});
    
    break_array r1, r2;
    if (breakpoints.front() != 0.0)
        {
            r1.emplace_back(std::make_pair(start, breakpoints.front()));
        }
    for (unsigned j = 1; j < breakpoints.size(); ++j)
        {
            double a = breakpoints[j - 1];
            double b = (j < breakpoints.size() - 1) ? breakpoints[j] : stop;
            if (j % 2 == 0.)
                {
                    r1.emplace_back(a, b);
                }
            else
                {
                    r2.emplace_back(a, b);
                }
        }
    return std::make_pair(std::move(r1), std::move(r2));
}

fwdpp::uint_t
ancestry_rec_mut_details(
    fwdpy11::SlocusPop& pop, ancestry_tracker& ancestry,
    std::queue<std::size_t>& gamete_recycling_bin,
    const fwdpp::uint_t parental_gamete1, const fwdpp::uint_t parental_gamete2,
    std::vector<double>& breakpoints,
    const std::vector<fwdpp::uint_t>& new_mutations,
    const std::tuple<ancestry_tracker::integer_type,
                     ancestry_tracker::integer_type>& pid,
    const ancestry_tracker::integer_type offspring_index)
{
	// This has the effect of removing any double x-overs.
	// Internally, fwdpp will do the right thing, but leaving
	// them in causes msprime to throw an error.  We could filter them
	// later, but doing it here is less code.
    breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()),
                      breakpoints.end());
    //split breakpoints among parental chromosomes and add edges
    auto breaks_pchrom = split_breakpoints(breakpoints);
    ancestry.add_edges(breaks_pchrom.first, std::get<0>(pid), offspring_index);
    ancestry.add_edges(breaks_pchrom.second, std::get<1>(pid), offspring_index);
    //add mutations
    ancestry.add_mutations(new_mutations, pop.mutations, offspring_index);
    return fwdpp::mutate_recombine(new_mutations, breakpoints, parental_gamete1, parental_gamete2,
                                     pop.gametes, pop.mutations, gamete_recycling_bin,
                                     pop.neutral, pop.selected);
}
