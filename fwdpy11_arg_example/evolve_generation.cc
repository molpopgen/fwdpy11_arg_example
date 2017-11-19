#include <fwdpp/diploid.hh>
#include "handle_recombination.hpp"

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
    ancestry_tracker& ancestry)
{

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
                pid, std::get<0>(offspring_indexes));
            breakpoints = recmodel(pop.gametes[p2g1], pop.gametes[p2g2],
                                   pop.mutations);
            pid = ancestry.get_parent_ids(p2, swap2);

            dip.second = ancestry_recombination_details(
                pop, ancestry, gamete_recycling_bin, p2g1, p2g2, breakpoints,
                pid, std::get<1>(offspring_indexes));

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
                                          pop.mcounts, 2 * N_next, std::true_type());
    // This is constant-time
    pop.diploids.swap(offspring);
}
