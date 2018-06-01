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
#include <fwdpy11/rng.hpp>
#include "ancestry_tracker.hpp"
#include <fwdpp/diploid.hh>
#include "handle_mut_rec.hpp"
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/fitness_models.hpp>
 
template <typename mcount_vec>
fwdpp::fwdpp_internal::recycling_bin_t<typename mcount_vec::size_type>
make_mut_queue(const mcount_vec &mcounts, ancestry_tracker& ancestry)
{
	fwdpp::fwdpp_internal::recycling_bin_t<typename mcount_vec::size_type> rv;
	const auto msize = mcounts.size();
	for (typename mcount_vec::size_type i = 0; i < msize; ++i)
		{
			if (!mcounts[i] && ancestry.preserve_mutation_index.find(i) == ancestry.preserve_mutation_index.end())
				rv.push(i);
		}
	return rv;
}

inline fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
w(fwdpy11::SlocusPop& pop)
{
    auto N_curr = pop.diploids.size();
    std::vector<double> fitnesses(N_curr);
    for (size_t i = 0; i < N_curr; ++i)
        {
            fitnesses[i] = pop.diploids[i].w; 
            pop.gametes[pop.diploids[i].first].n = 0;
            pop.gametes[pop.diploids[i].second].n = 0;
        }
    auto lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
        gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
    return lookup;
}

template <typename breakpoint_function, typename mutation_model>
void
evolve_generation(
    const fwdpy11::GSLrng_t& rng, fwdpy11::SlocusPop & pop,
    const fwdpp::uint_t N_next, const double mu,
     const mutation_model& mmodel, const breakpoint_function & recmodel,
     ancestry_tracker& ancestry)
{

    auto gamete_recycling_bin
        = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
    auto mutation_recycling_bin
        = make_mut_queue(pop.mcounts, ancestry);

    // Efficiency hit.  Unavoidable
    // in use case of a sampler looking
    // at the gametes themselves (even tho
    // gamete.n has little bearing on anything
    // beyond recycling).  Can revisit later
    for (auto&& g : pop.gametes)
        g.n = 0;

	auto lookup = w(pop);
	auto ff = fwdpp::multiplicative_diploid();
    decltype(pop.diploids) offspring(N_next);

    // Generate the offspring
    std::size_t label = 0;
    for (auto& dip : offspring)
        {
			auto p1 = gsl_ran_discrete(rng.get(), lookup.get());
            auto p2 = gsl_ran_discrete(rng.get(), lookup.get());
            auto p1g1 = pop.diploids[p1].first;
            auto p1g2 = pop.diploids[p1].second;
            auto p2g1 = pop.diploids[p2].first;
            auto p2g2 = pop.diploids[p2].second;

            // Mendel
            int swap1 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
            int swap2 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
            if (swap1)
                std::swap(p1g1, p1g2);
            if (swap2)
                std::swap(p2g1, p2g2);

			auto breakpoints = recmodel();
			
            auto new_mutations = fwdpp::generate_new_mutations(
                mutation_recycling_bin, rng.get(), mu, pop.diploids[p1], pop.gametes, pop.mutations, p1g1, mmodel);
            auto pid = ancestry.get_parent_ids(p1, swap1);
            auto offspring_indexes = ancestry.get_next_indexes(0);
            dip.first = ancestry_rec_mut_details(
                pop, ancestry, gamete_recycling_bin, p1g1, p1g2, breakpoints, new_mutations,
                pid, std::get<0>(offspring_indexes));
                
			breakpoints = recmodel();
            new_mutations = fwdpp::generate_new_mutations(
                mutation_recycling_bin, rng.get(), mu, pop.diploids[p1], pop.gametes, pop.mutations, p2g1, mmodel);
            pid = ancestry.get_parent_ids(p2, swap2);
            dip.second = ancestry_rec_mut_details(
                pop, ancestry, gamete_recycling_bin, p2g1, p2g2, breakpoints, new_mutations,
                pid, std::get<1>(offspring_indexes));

            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;

            assert(pop.gametes[dip.first].n);
            assert(pop.gametes[dip.second].n);
            dip.label = label++;
            dip.deme = 0;
            dip.sex = 0;
            dip.parental_data = std::make_tuple(p1,p2);
            dip.w = ff(dip, pop.gametes, pop.mutations);
        }
    ancestry.finish_generation();
    fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                           pop.mcounts);
    fwdpp::fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations,
                                          pop.mcounts, 2 * N_next, std::true_type());
    // This is constant-time
    pop.diploids.swap(offspring);
}



#endif
