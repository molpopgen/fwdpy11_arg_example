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
#include <gsl/gsl_randist.h>
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

struct parent_lookup_tables
// Our object for choosing parents each generation
{
    // These return indexes of parents from demes 1 and 2,
    // resp, chosen in O(1) time proportional to
    // relative fitness within each deme
    fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr lookup1, lookup2;
    // These vectors map indexes returned from sampling
    // lookup1 and lookup2 to diploids in the population
    // object.
    std::vector<std::size_t> parents1, parents2;
};

inline parent_lookup_tables
migrate_and_calc_fitness(const gsl_rng *r, fwdpy11::SlocusPop& pop, 
	                     const fwdpp::uint_t N1, const fwdpp::uint_t N2, 
	                     const double m12, const double m21, const bool expect_split)
// This function will be called at the start of each generation.
// The main goal is to return the lookup tables described above.
// But, "while we're at it", it does some other stuff that
// needs to be done at the start of each generation.
// Neither the most rigorous nor the most efficient:
// 1. Ignores probability of back-migration.
// 2. Allocates 4 vectors each generation.
{
    parent_lookup_tables rv;

    // Temp containers for fitnesses in each deme,
    // post-migration
    std::vector<double> w1, w2;

    // Pick no. migrants 1 -> 2 and 2 -> 1.
    unsigned nmig12 = gsl_ran_poisson(r, static_cast<double>(N1) * m12);
    unsigned nmig21 = gsl_ran_poisson(r, static_cast<double>(N2) * m21);

	if(N2 == 0 && m12 > 0){ nmig12 = expect_split ? static_cast<double>(N1) * m12 : std::max(nmig12,2U); } //population 2 is being established, must have at least 2 individuals migrating to it
	nmig12 = std::min(nmig12,N1); //can't have more individuals migrating than present in each population
	nmig21 = std::min(nmig21,N2);
	
    // Fill a vector of N1 zeros and N2 ones:
    std::vector<fwdpp::uint_t> deme_labels(N1, 0);
    deme_labels.resize(N1 + N2, 1);
    assert(deme_labels.size() == pop.diploids.size());

    // Set up source and destination containers
    // for sampling w/o replacement
    std::vector<std::size_t> individuals(N1 + N2);
    std::iota(std::begin(individuals), std::end(individuals), 0);
    std::vector<std::size_t> migrants(std::max(nmig12, nmig21));

    // Who is migrating 1 -> 2?
    gsl_ran_choose(r, migrants.data(), nmig12, individuals.data(), N1,
                   sizeof(std::size_t));
    for (std::size_t i = 0; i < nmig12; ++i)
        {
            deme_labels[migrants[i]] = !deme_labels[migrants[i]];
        }

    // Exact same logic for migrants 2 -> 1
    gsl_ran_choose(r, migrants.data(), nmig21, individuals.data() + N1, N2,
                   sizeof(std::size_t));
    for (std::size_t i = 0; i < nmig21; ++i)
        {
            deme_labels[migrants[i]] = !deme_labels[migrants[i]];
        }

    // Go over all parents, set gametes counts to zero,
    // and put individual IDs and fitnesses into
    // the right vectors:
    for (std::size_t i = 0; i < deme_labels.size(); ++i)
        {
            // fwdpp requires that we zero out gamete
            // counts each generation.  Since we're looping
            // over diploids here, now is a good time to
            // handle this task, which saves us from having to
            // do another O(N1+N2) loop:
            pop.gametes[pop.diploids[i].first].n
                = pop.gametes[pop.diploids[i].second].n = 0;
            if (deme_labels[i] == 0)
                {
                    rv.parents1.push_back(i);
                    w1.push_back(pop.diploid_metadata[i].w);
                }
            else
                {
                    rv.parents2.push_back(i);
                    w2.push_back(pop.diploid_metadata[i].w);
                }
        }

    // Set up our lookup tables:
    rv.lookup1.reset(gsl_ran_discrete_preproc(rv.parents1.size(), w1.data()));
    rv.lookup2.reset(gsl_ran_discrete_preproc(rv.parents2.size(), w2.data()));
    return rv;
};

inline parent_lookup_tables
duplicate_and_calc_fitness(const gsl_rng *r, fwdpy11::SlocusPop& pop, 
	                     const fwdpp::uint_t N1)
{
    parent_lookup_tables rv;
    // Temp containers for fitnesses in each deme,
    // post-migration
    std::vector<double> w1, w2;

    // Go over all parents, set gametes counts to zero,
    // and put individual IDs and fitnesses into
    // the right vectors:
    for (std::size_t i = 0; i < N1; ++i)
        {
            // fwdpp requires that we zero out gamete
            // counts each generation.  Since we're looping
            // over diploids here, now is a good time to
            // handle this task, which saves us from having to
            // do another O(N1+N2) loop:
            pop.gametes[pop.diploids[i].first].n
                = pop.gametes[pop.diploids[i].second].n = 0;
            rv.parents1.push_back(i);
            w1.push_back(pop.diploid_metadata[i].w);
            rv.parents2.push_back(i);
            w2.push_back(pop.diploid_metadata[i].w);
        }

    // Set up our lookup tables:
    rv.lookup1.reset(gsl_ran_discrete_preproc(rv.parents1.size(), w1.data()));
    rv.lookup2.reset(gsl_ran_discrete_preproc(rv.parents2.size(), w2.data()));
    return rv;
};

inline parent_lookup_tables
migration_fitness_parents(const gsl_rng *r, fwdpy11::SlocusPop& pop, 
	                     const fwdpp::uint_t prev_N1, const fwdpp::uint_t prev_N2, 
	                     const double m12, const double m21, const bool expect_split)
{
	if(m12 == 1.f && prev_N2 == 0)
		{ 
			return duplicate_and_calc_fitness(r, pop, prev_N1); 
	   	}
	else
	   	{
	    	return migrate_and_calc_fitness(r, pop, prev_N1, prev_N2, m12, m21, expect_split);
	   	}

}

inline auto 
get_parent(const gsl_rng *r, const parent_lookup_tables& lookups, int deme) -> decltype(lookups.parents1[0])
{
	if (deme == 0) // pick parent from pop 1
		{
			return lookups.parents1[gsl_ran_discrete(r, lookups.lookup1.get())];
		}
	else // pick parent from pop 2
		{
			return lookups.parents2[gsl_ran_discrete(r, lookups.lookup2.get())];
		}
}

template <typename fitness_fxn, typename mutation_fxn, typename breakpoint_fxn>
void
evolve_generation(
    const fwdpy11::GSLrng_t& rng, fwdpy11::SlocusPop& pop,
    const fwdpp::uint_t N1, const fwdpp::uint_t prev_N2, const fwdpp::uint_t N2, 
    const double m12, const double m21, fwdpp::uint_t & split_N1_loss, const bool recover_split_loss, const bool expect_split, const double mu, const fitness_fxn& wmodel, 
    const mutation_fxn& mmodel, const breakpoint_fxn& rmodel, 
    ancestry_tracker& ancestry)
{

    auto gamete_recycling_bin
        = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
    auto mutation_recycling_bin
        = make_mut_queue(pop.mcounts, ancestry);
    auto lookups 
		= migration_fitness_parents(rng.get(), pop, pop.N, prev_N2, m12, m21, expect_split);
		
	if(m12 < 1.f && prev_N2 == 0 && !recover_split_loss){
		split_N1_loss = pop.N - lookups.parents1.size(); //allows split size to be random and decrease N1 by exact amount of individuals who left the population to form N2
	}
	
	auto N1_loss = N1 - split_N1_loss; 
	
    decltype(pop.diploids) offspring(N1_loss+N2);
    decltype(pop.diploid_metadata) offspring_metadata(N1_loss+N2);
    
    // Generate the offspring
    std::size_t label = 0;
    for (auto& dip : offspring)
        {
            int deme = (label >= N1_loss);
            auto offspring_indexes = ancestry.get_next_indexes(deme);
        	
	    	auto p1 = get_parent(rng.get(), lookups, deme);
            auto p2 = get_parent(rng.get(), lookups, deme);
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
			
	   		auto breakpoints = rmodel();
            auto new_mutations = fwdpp::generate_new_mutations(
                mutation_recycling_bin, rng.get(), mu, pop.diploids[p1], pop.gametes, pop.mutations, p1g1, mmodel);
            auto pid = ancestry.get_parent_ids(p1, swap1);
            dip.first = ancestry_rec_mut_details(
                pop, ancestry, gamete_recycling_bin, p1g1, p1g2, breakpoints, new_mutations,
                pid, std::get<0>(offspring_indexes));
                
			breakpoints = rmodel();
            new_mutations = fwdpp::generate_new_mutations(
                mutation_recycling_bin, rng.get(), mu, pop.diploids[p1], pop.gametes, pop.mutations, p2g1, mmodel);
            pid = ancestry.get_parent_ids(p2, swap2);
            dip.second = ancestry_rec_mut_details(
                pop, ancestry, gamete_recycling_bin, p2g1, p2g2, breakpoints, new_mutations,
                pid, std::get<1>(offspring_indexes));

            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;

            offspring_metadata[label].label = label;
            offspring_metadata[label].deme = deme;
            offspring_metadata[label].sex = 0;
            offspring_metadata[label].parents[0] = p1;
            offspring_metadata[label].parents[1] = p2;
            offspring_metadata[label].w = wmodel(dip, pop.gametes, pop.mutations);
            label++;
        }
    ancestry.finish_generation();
    fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                           pop.mcounts);
    fwdpp::fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations,
                                          pop.mcounts, 2 * (N1_loss+N2), std::true_type());
    // This is constant-time
    pop.diploids.swap(offspring);
    pop.diploid_metadata.swap(offspring_metadata);
}



#endif
