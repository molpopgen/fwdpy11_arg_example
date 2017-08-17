#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <fwdpy11/rules/wf_rules.hpp>
#include <fwdpy11/sim_functions.hpp>
#include "evolve_generation.hpp"
#include "ancestry_tracker.hpp"

namespace py = pybind11;

void
evolve_singlepop_regions_track_ancestry(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
	py::function ancestry_processor,
    py::array_t<std::uint32_t> popsizes, const double mu_selected,
    const double recrate, const KTfwd::extensions::discrete_mut_model& mmodel,
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
    auto rules = fwdpy11::wf_rules();

    auto fitness_callback = fitness.callback();
    fitness.update(pop);
    ancestry_tracker ancestry(pop.N);
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
			py::print("generation = ",generation);
			py::bool_ processor_rv = ancestry_processor(pop.generation,ancestry.nodes,ancestry.edges);
			bool gc = processor_rv.cast<bool>();
			py::print("GC = ", gc);
            const auto N_next = popsizes.at(generation);
            evolve_generation(
                rng, pop, N_next, mu_selected, mmodels, recmap,
                std::bind(&fwdpy11::wf_rules::pick1, &rules,
                          std::placeholders::_1, std::placeholders::_2),
                std::bind(&fwdpy11::wf_rules::pick2, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, selfing_rate),
                std::bind(&fwdpy11::wf_rules::update, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, std::placeholders::_4,
                          std::placeholders::_5),
                ancestry, std::true_type());
			py::print("evolved!");
            pop.N = N_next;
            fwdpy11::update_mutations_wrapper()(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N);
            fitness.update(pop);
        }
    --pop.generation;
}
