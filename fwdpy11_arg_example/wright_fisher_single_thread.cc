#include <chrono>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/chrono.h>
#include <fwdpy11/rules/wf_rules.hpp>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpp/extensions/regions.hpp>
#include "evolve_generation.hpp"

namespace py = pybind11;

template <typename mcont_t, typename fixation_container_t,
		  typename fixation_time_container_t,
		  typename mutation_lookup_table>
void
update_mutations(mcont_t &mutations, fixation_container_t &fixations,
				 fixation_time_container_t &fixation_times,
				 mutation_lookup_table &lookup, ancestry_tracker & ancestry,
				 std::vector<KTfwd::uint_t> &mcounts,
				 const unsigned &generation, const unsigned &twoN,
				 const bool remove_selected_fixations);

// This function runs the simulation itself.
// The details of a generation are in the file
// evolve_generation.hpp. This function gets
// exposed to Python, but it is not called
// directly by a user. The function "evolve_track"
// in evolve_arg.py takes the fwdpy11 objects and
// passes them along to this function. The main purpose
// of this function is to check parameters and send
// data to msprime when needed.
// The argument ancestry_processor is a Python callable.
// It is to handle the GC/simplification step via msprime.
// It should be an instance of ARGsimplifier.
// The return value is the time spent simulating.
double
evolve_singlepop_regions_track_ancestry(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,  ancestry_tracker & ancestry, 
    py::function arg_simplifier, pybind11::function anc_sampler, py::array_t<std::uint32_t> popsizes, 
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
        
    const auto generations = sample_simplify.attr("total_generations").cast<unsigned>();
    
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
    anc_sampler();
    ++pop.generation;
    auto rules = fwdpy11::wf_rules();

    auto fitness_callback = fitness.callback();
    fitness.update(pop);
    auto wbar = rules.w(pop, fitness_callback);

    double time_simulating = 0.0;
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            const auto N_next = popsizes.at(generation);
            auto start = std::clock();
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
                ancestry);
            pop.N = N_next;
            update_mutations(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, ancestry, pop.mcounts, pop.generation, 2 * pop.N, false);
            fitness.update(pop);
            wbar = rules.w(pop, fitness_callback);
            auto stop = std::clock();
            auto dur = (stop - start) / static_cast<double>(CLOCKS_PER_SEC);
            time_simulating += dur;
            //Ask if we need to garbage collect:
            arg_simplifier(anc_sampler());
        }
        
    arg_simplifier(true);
    --pop.generation;
    return time_simulating;
}

template <typename mcont_t, typename fixation_container_t,
		  typename fixation_time_container_t,
		  typename mutation_lookup_table>
void
update_mutations(mcont_t &mutations, fixation_container_t &fixations,
				 fixation_time_container_t &fixation_times,
				 mutation_lookup_table &lookup, ancestry_tracker & ancestry,
				 std::vector<KTfwd::uint_t> &mcounts,
				 const unsigned &generation, const unsigned &twoN,
				 const bool remove_selected_fixations)
{
	using namespace KTfwd;
	static_assert(
		typename traits::is_mutation_t<
			typename mcont_t::value_type>::type(),
		"mutation_type must be derived from KTfwd::mutation_base");
	assert(mcounts.size() == mutations.size());
	for (unsigned i = 0; i < mcounts.size(); ++i)
		{
			assert(mcounts[i] <= twoN);
			if (mcounts[i] == twoN)
				{
					auto loc = std::lower_bound(
						fixations.begin(), fixations.end(),
						mutations[i].pos,
						[](const typename fixation_container_t::value_type
							   &__mut,
						   const double &__value) noexcept {
							return __mut.pos < __value;
						});
					auto d = std::distance(fixations.begin(), loc);
					if (mutations[i].neutral
						|| remove_selected_fixations == true)
						{
							fixations.insert(loc, mutations[i]);
							fixation_times.insert(
								fixation_times.begin() + d, generation);
							mcounts[i] = 0; // set count to zero to mark
											// mutation as "recyclable"
							lookup.erase(mutations[i].pos); // remove
															// mutation
															// position
															// from lookup
						}
					else
						{
							if (loc == fixations.end()
								|| (loc->pos != mutations[i].pos
									&& loc->g != mutations[i].g))
								{
									fixations.insert(loc, mutations[i]);
									fixation_times.insert(
										fixation_times.begin() + d,
										generation);
								}
						}
				}
			if (!mcounts[i] && ancestry.preserve_mutation_index.find(i) == ancestry.preserve_mutation_index.end())
				lookup.erase(mutations[i].pos);
		}
}
