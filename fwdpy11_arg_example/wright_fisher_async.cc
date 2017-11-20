#include <chrono>
#include <future>
#include <pybind11/chrono.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <fwdpy11/rules/wf_rules.hpp>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpp/extensions/regions.hpp>
#include "evolve_generation.hpp"
#include "wf_rules_async.hpp"

namespace py = pybind11;

double
evolve_singlepop_regions_track_ancestry_async(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    ancestry_tracker& ancestry, py::function ancestry_processor,
    const KTfwd::uint_t gc_interval, py::array_t<std::uint32_t> popsizes,
    const double mu_selected, const double recrate,
    const KTfwd::extensions::discrete_mut_model& mmodel,
    const KTfwd::extensions::discrete_rec_model& rmodel,
    fwdpy11::single_locus_fitness& fitness, const double selfing_rate)
{
    if (pop.generation > 0)
        {
            throw std::runtime_error(
                "this population has already been evolved.");
        }
    const auto generations = popsizes.size();
    if (!generations)
        throw std::runtime_error("empty list of population sizes");
    if (mu_selected < 0.)
        {
            throw std::runtime_error("negative selected mutation rate: "
                                     + std::to_string(mu_selected));
        }
    if (recrate < 0.)
        {
            throw std::runtime_error("negative recombination rate: "
                                     + std::to_string(recrate));
        }
    pop.mutations.reserve(
        std::ceil(std::log(2 * pop.N)
                  * (4. * double(pop.N) * (mu_selected)
                     + 0.667 * (4. * double(pop.N) * (mu_selected)))));
    const auto recmap = KTfwd::extensions::bind_drm(
        rmodel, pop.gametes, pop.mutations, rng.get(), recrate);
    const auto mmodels = KTfwd::extensions::bind_dmm(
        mmodel, pop.mutations, pop.mut_lookup, rng.get(), 0.0, mu_selected,
        &pop.generation);
    ++pop.generation;
    auto rules = wf_rules_async_fitness(2);

    auto fitness_callback = fitness.callback();
    fitness.update(pop);
    auto wbar = rules.w(pop, fitness_callback);

    std::future<py::object> msprime_future;
    ancestry_tracker local_ancestry_tracker;
    double time_simulating = 0.0;
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            if (pop.generation > 0 && pop.generation % gc_interval == 0.)
                {
                    if (msprime_future.valid())
                        {
                            msprime_future.wait();
                            auto result = msprime_future.get();
                            auto result_tuple = result.cast<py::tuple>();
                            ancestry.post_process_gc(result_tuple, false);
                        }

                    ancestry.exchange_for_async(local_ancestry_tracker);
                    msprime_future = std::async(
                        std::launch::async, ancestry_processor, pop.generation,
                        std::ref(local_ancestry_tracker));
                }
            //This is not great API design, but
            //we need to clear the offspring indexes here:
            ancestry.offspring_indexes.clear();
            const auto N_next = popsizes.at(generation);
            auto start = std::clock();
            evolve_generation(
                rng, pop, N_next, mu_selected, mmodels, recmap,
                std::bind(&decltype(rules)::pick1, &rules,
                          std::placeholders::_1, std::placeholders::_2),
                std::bind(&decltype(rules)::wf_rules::pick2, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, selfing_rate),
                std::bind(&decltype(rules)::update, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, std::placeholders::_4,
                          std::placeholders::_5),
                ancestry);
            pop.N = N_next;
            fwdpy11::update_mutations(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N, true);
            fitness.update(pop);
            wbar = rules.w(pop, fitness_callback);
            auto stop = std::clock();
            auto dur = (stop - start) / static_cast<double>(CLOCKS_PER_SEC);
            time_simulating += dur;
        }
    --pop.generation;
    if (msprime_future.valid())
        {
            msprime_future.wait();
            auto result = msprime_future.get();
            auto result_tuple = result.cast<py::tuple>();
            ancestry.post_process_gc(result_tuple, false);
        }

    return time_simulating;
}
