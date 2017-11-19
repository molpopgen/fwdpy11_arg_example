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
evolve_singlepop_regions_track_ancestry_python_queue(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    ancestry_tracker& ancestry, py::object python_queue,
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

    double time_simulating = 0.0;

    std::size_t python_qsize
        = python_queue.attr("maxsize").cast<std::size_t>();
    py::gil_scoped_release GIL_release;
    std::vector<py::object> faux_memory_pool(python_qsize);
    for (auto& i : faux_memory_pool)
        {
            i = py::cast(ancestry_tracker(ancestry));
        }
    std::size_t items_submitted = 0;
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            if (pop.generation > 0 && pop.generation % gc_interval == 0.)
                {
                    {
                        py::gil_scoped_acquire acquire;
                        ancestry.exchange_for_async(
                            faux_memory_pool[items_submitted]
                                .cast<ancestry_tracker&>());
                        python_queue.attr("put")(py::make_tuple(
                            pop.generation,
                            faux_memory_pool[items_submitted++]));
                    }
                    if (items_submitted >= python_qsize)
                        {
                            items_submitted = 0;
                        }
                    ancestry.first_parental_index = 0;
                    ancestry.next_index = 2 * pop.diploids.size();
                    ancestry.nodes.clear();
                    ancestry.edges.clear();
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

    return time_simulating;
}
