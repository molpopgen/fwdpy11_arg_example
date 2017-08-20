#ifndef FWDPY11_ARG_EXAMPLE_EVOLVE_GENERATION_HPP__
#define FWDPY11_ARG_EXAMPLE_EVOLVE_GENERATION_HPP__

#include <cstddef>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <queue>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <fwdpp/diploid.hh>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <fwdpy11/rng.hpp>
#include "ancestry_edge_sets.hpp"
#include "handle_recombination.hpp"

// This is a copy/paste + modification of fwdpy11's
// existing function to evolve a single-deme, single-region
// population object for one generation.  The modification is to include updating
// an ancestry_tracker object
template <typename poptype, typename pick1_function, typename pick2_function,
          typename update_function, typename mutation_model,
          typename recombination_model, typename mutation_removal_policy>
void
evolve_generation(const fwdpy11::GSLrng_t& rng, poptype& pop,
                  const KTfwd::uint_t N_next, const double mu,
                  const mutation_model& mmodel,
                  const recombination_model& recmodel,
                  const pick1_function& pick1, const pick2_function& pick2,
                  const update_function& update, ancestry_tracker& ancestry,
                  const mutation_removal_policy& mrp)
{
    static_assert(std::is_same<typename poptype::popmodel_t,
                               KTfwd::sugar::SINGLEPOP_TAG>::value,
                  "Population type must be a single-locus, single-deme type.");

    auto gamete_recycling_bin
        = KTfwd::fwdpp_internal::make_gamete_queue(pop.gametes);
    auto mutation_recycling_bin
        = KTfwd::fwdpp_internal::make_mut_queue(pop.mcounts);

    // Efficiency hit.  Unavoidable
    // in use case of a sampler looking
    // at the gametes themselves (even tho
    // gamete.n has little bearing on anything
    // beyond recycling).  Can revisit later
    for (auto&& g : pop.gametes)
        g.n = 0;

    decltype(pop.diploids) offspring(N_next);

    // Generate the offspring
    std::size_t label = 0;
    for (auto& dip : offspring)
        {
            auto p1 = pick1(rng, pop);
            auto p2 = pick2(rng, pop, p1);
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

            auto breakpoints = recmodel(pop.gametes[p1g1], pop.gametes[p1g2],
                                        pop.mutations);
            auto pid = ancestry.get_parent_ids(p1, swap1);
            auto offspring_indexes = ancestry.get_next_indexes();
            dip.first = ancestry_recombination_details(
                pop, ancestry, gamete_recycling_bin, p1g1, p1g2, breakpoints,
                pid, offspring_indexes);
            breakpoints = recmodel(pop.gametes[p2g1], pop.gametes[p2g2],
                                   pop.mutations);
            pid = ancestry.get_parent_ids(p2, swap2);
            offspring_indexes = ancestry.get_next_indexes();
            dip.second = ancestry_recombination_details(
                pop, ancestry, gamete_recycling_bin, p2g1, p2g2, breakpoints,
                pid, offspring_indexes);

            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;

            // now, add new mutations
            dip.first = KTfwd::mutate_gamete_recycle(
                mutation_recycling_bin, gamete_recycling_bin, rng.get(), mu,
                pop.gametes, pop.mutations, dip.first, mmodel,
                KTfwd::emplace_back());
            dip.second = KTfwd::mutate_gamete_recycle(
                mutation_recycling_bin, gamete_recycling_bin, rng.get(), mu,
                pop.gametes, pop.mutations, dip.second, mmodel,
                KTfwd::emplace_back());

            assert(pop.gametes[dip.first].n);
            assert(pop.gametes[dip.second].n);
            dip.label = label++;
            update(rng, dip, pop, p1, p2);
        }
    ancestry.finish_generation();
    KTfwd::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                           pop.mcounts);
    KTfwd::fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations,
                                          pop.mcounts, 2 * N_next, mrp);
    // This is constant-time
    pop.diploids.swap(offspring);
}

#endif
