#include <chrono>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/chrono.h>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpy11/policies/mutation.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>
#include "evolve_generation.hpp"

namespace py = pybind11;

template <typename mcont_t, typename fixation_container_t,
		  typename fixation_time_container_t,
		  typename mutation_lookup_table>
void
update_mutations(mcont_t &mutations, fixation_container_t &fixations,
				 fixation_time_container_t &fixation_times,
				 mutation_lookup_table &lookup, ancestry_tracker & ancestry,
				 std::vector<fwdpp::uint_t> &mcounts,
				 const unsigned &generation, const unsigned &twoN);

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
evolve_track_ancestry(
    const fwdpy11::GSLrng_t& rng, fwdpy11::SlocusPop & pop,  ancestry_tracker & ancestry, 
    py::function sample_simplify, py::array_t<std::uint32_t> popsizes, 
    py::array_t<std::uint32_t> pop2array, py::array_t<float> migarray, 
    const double mu_selected, 
    const fwdpp::extensions::discrete_mut_model<fwdpy11::SlocusPop::mcont_t>
        &mmodel,
    const fwdpp::extensions::discrete_rec_model &rmodel)
{
    if (pop.generation > 0)
        {
            throw std::runtime_error(
                "this population has already been evolved.");
        }
        
    const auto generations = sample_simplify.attr("total_generations").cast<unsigned>();
    
    if (!generations)
        throw std::runtime_error("empty list of population sizes");
    if(!pop2array)
    	throw std::runtime_error("empty population 2 array");
    if(!migarray)
    	throw std::runtime_error("empty migration array");
    if (mu_selected < 0.)
        {
            throw std::runtime_error("negative selected mutation rate: "
                                     + std::to_string(mu_selected));
        }
    
    pop.mutations.reserve(
        std::ceil(std::log(2 * pop.N)
                  * (4. * double(pop.N) * (mu_selected)
                     + 0.667 * (4. * double(pop.N) * (mu_selected)))));
        
    const auto bound_mmodel = fwdpp::extensions::bind_dmm(rng.get(), mmodel);
    const auto bound_rmodel = [&rng, &rmodel]() { return rmodel(rng.get()); };
    
    const auto ff = fwdpp::multiplicative_diploid(fwdpp::fitness(2.0));   
    
    ++pop.generation; //simulation starts at demography generation 1
    auto prev_pop2size = 0U;
	const auto pop2_start = pop2array.at(1);
	const auto pop2_end = pop2array.at(2);
	const auto mig_start = migarray.at(3);
	const auto mig_end = migarray.at(4);
	const bool recover_split = migarray.at(5);
    double time_simulating = 0.0;
    
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
        	auto pop2size = 0U;
        	auto mig12 = 0.f;
        	auto mig21 = 0.f;
        	if(pop.generation >= pop2_start && pop.generation < pop2_end){
        		pop2size = pop2array.at(0);
        	}
        	if(pop.generation >= mig_start && pop.generation < mig_end){
        		mig12 = migarray.at(0);
        		mig21 = migarray.at(1);
        	}
        	if(pop.generation == pop2_start){
        		mig12 = migarray.at(2);
        		mig21 = 0;
        	}
            const auto N_next = popsizes.at(pop.generation); //popsizes.at(0) is pop.N
            auto start = std::clock();
            evolve_generation(
                rng, pop, N_next, prev_pop2size, pop2size, mig12, mig21, recover_split, mu_selected, ff, bound_mmodel, bound_rmodel, ancestry);
            pop.N = N_next;
            prev_pop2size = pop2size;
            update_mutations(
                pop.mutations, pop.fixations, pop.fixation_times, pop.mut_lookup, 
                ancestry, pop.mcounts, pop.generation, 2 * (pop.N+prev_pop2size));
            auto stop = std::clock();
            auto dur = (stop - start) / static_cast<double>(CLOCKS_PER_SEC);
            time_simulating += dur;
            //Ask if we need to sample and/or garbage collect:
            sample_simplify();
        }
        
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
				 std::vector<fwdpp::uint_t> &mcounts,
				 const unsigned &generation, const unsigned &twoN)
{
	using namespace fwdpp;
	static_assert(
		typename traits::is_mutation_t<
			typename mcont_t::value_type>::type(),
		"mutation_type must be derived from fwdpp::mutation_base");
	assert(mcounts.size() == mutations.size());
	for (unsigned i = 0; i < mcounts.size(); ++i)
		{
			assert(mcounts[i] <= twoN);
			if (mcounts[i] == twoN)
				{
					fixations.push_back(mutations[i]);
					fixation_times.push_back(generation);
					ancestry.preserve_fixation(i);
				}
			else if (!mcounts[i] && ancestry.preserve_mutation_index.find(i) == ancestry.preserve_mutation_index.end())
				{
					auto itr = lookup.equal_range(mutations[i].pos);
					if (itr.first != lookup.end())
						{
							while (itr.first != itr.second)
								{
									if (itr.first->second == i)
										{
											lookup.erase(itr.first);
											break;
										}
									++itr.first;
								}
						}
				}
		}
}
