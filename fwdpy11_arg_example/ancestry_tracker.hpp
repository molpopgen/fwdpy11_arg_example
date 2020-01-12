// This file defines a C++ class to record
// nodes and edges during a sim.  It also
// gets the data ready to send to msprime
// and cleanup after msprime simplifies things.
// The ancestry_tracer is mainly a book-keeper
// and it is exposed to Python as
// wfarg.AncestryTracker.
//
// This implementation does not have to be
// header-only.  However, it is easier (lazier)
// to do that, and we'll fix that when moving things to
// fwdpy11.

#ifndef ANCESTRY_ANCESTRY_TRACKER_HPP__
#define ANCESTRY_ANCESTRY_TRACKER_HPP__

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <limits>
#include <cstdint>
#include <atomic>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fwdpy11/types/SlocusPop.hpp>

#include "node.hpp"
#include "edge.hpp"
#include "mutation.hpp"

// struct ancestry_data
// {
//     using integer_type = decltype(edge::parent);
//     /// Nodes:
//     std::vector<node> nodes;
//     /// The ARG:
//     std::vector<edge> edges;
//     /// Mutations:
//     std::vector<mutation> mutations;
//     std::vector<integer_type> samples;
//     pybind11::object lock_;
//     ancestry_data() : nodes{}, edges{}, mutations{}, samples{}, lock_{}
//     {
//         pybind11::module threading = pybind11::module::import("threading");
//         lock_ = threading.attr("Lock")();
//     }
// };

struct ancestry_tracker
{
    using integer_type = decltype(edge::parent);
    using index_pair = std::pair<integer_type,integer_type>;
    using mut_index_set = std::unordered_set<std::uint32_t>;
    /// Nodes:
    std::vector<node> nodes;
    /// The ARG:
    std::vector<edge> edges;
    /// The edges generated for each generation:
    std::vector<edge> temp;
    /// Mutations:
    std::vector<mutation> mutations;
    /// start-end node IDs for a generation:
    index_pair node_indexes;
    /// indices of mutations to preserve in the simulation
    mut_index_set preserve_mutation_index;
    /// current generation, total generation, next node ID to use, current index generation
    integer_type generation, total_generations, next_index;
    ancestry_tracker(const integer_type N, 
                     const integer_type next_index_,
                     const integer_type total_generations_)
        : nodes{ }, edges{ }, temp{ }, mutations{ }, preserve_mutation_index{ }, 
          generation{ 1 }, total_generations{ total_generations_ },
          next_index{ next_index_ }
    {
        nodes.reserve(2 * N);
        edges.reserve(2 * N);
        temp.reserve(N);
        //no need to reserve mutation space
		
        //if next_index == 0, then did not initialize with tree sequence
        //so emplace back 2N nodes for the generation 0 and set next_index to 2N
        if (next_index == 0)
            {
                for (integer_type i = 0; i < 2 * N; ++i)
                    {
                        //ID, time 0, population 0
                        double rev_gen = total_generations;
                        nodes.emplace_back(node{i, 0, rev_gen});
                    }
                next_index = 2 * N;
            }
        node_indexes = std::make_pair(0,next_index);
    }

    std::tuple<integer_type, integer_type>
    get_parent_ids(const std::uint32_t p, const int did_swap)
    {
        return std::make_tuple(
            node_indexes.first + 2 * static_cast<integer_type>(p) + did_swap, // node_indexes.first == first_parental_index
            node_indexes.first + 2 * static_cast<integer_type>(p)+ !did_swap);
    }

    std::tuple<integer_type, integer_type>
    get_next_indexes(const std::int32_t population)
    {
        auto rv = std::make_tuple(next_index, next_index + 1);
        double rev_gen = total_generations - generation;
        nodes.emplace_back(node{std::get<0>(rv), population, rev_gen});
        nodes.emplace_back(node{std::get<1>(rv), population, rev_gen});
        next_index += 2;
        return rv;
    }

    void
    add_edges(const std::vector<std::pair<double, double>>& breakpoints,
              const integer_type parent, const integer_type child)
    {
        for (auto&& bi : breakpoints)
            {
                temp.emplace_back(edge{bi.first, bi.second, parent, child});
            }
    }

    template <typename mcont_t>
    void
    add_mutations(const std::vector<std::uint32_t>& mutation_ids, mcont_t &popmut, const integer_type node_id)
    {
        for (auto&& mut_id : mutation_ids)
            {
                mutations.emplace_back(mutation{node_id, popmut[mut_id].pos, mut_id});
            }
    }
    
    void
    preserve_mutations_sample(const pybind11::array_t<integer_type> & indiv_samples, fwdpy11::SlocusPop& pop)
    {
    	for(auto i : indiv_samples){
    		int index = i.cast<integer_type>();
    		auto g1 = pop.diploids[index].first;
    		auto g2 = pop.diploids[index].second;
    		
    		preserve_mutation_index.insert(pop.gametes[g1].smutations.begin(),pop.gametes[g1].smutations.end());
    		preserve_mutation_index.insert(pop.gametes[g2].smutations.begin(),pop.gametes[g2].smutations.end());
    	}
    }
    
    void
    preserve_fixation(const std::uint32_t & mutation_id)
    {
    	preserve_mutation_index.insert(mutation_id);
    }

    void
    finish_generation()
    {
        edges.insert(edges.end(), temp.begin(), temp.end());
        temp.clear();
        node_indexes.first = node_indexes.second;
        node_indexes.second = next_index; 
        ++generation;
    }
    
    void
    pre_process_gc(fwdpy11::SlocusPop& pop)
    {
        mutations.erase(std::remove_if(mutations.begin(),mutations.end(),[&pop](const mutation & m) { return (pop.mutations[m.mutation_id].pos != m.pos); }), mutations.end());
    }

    void
    post_process_gc(integer_type _next_index)
    {
        next_index = _next_index;
        node_indexes.first = 0;
        node_indexes.second = next_index; 
           
        nodes.clear();
        edges.clear();
        mutations.clear();
    }

//     void
//     exchange_for_async(ancestry_data& a)
//     {
//         nodes.swap(a.nodes);
//         edges.swap(a.edges);
//         mutations.swap(a.mutations);
//         a.samples.assign(offspring_indexes.begin(), offspring_indexes.end());
//         nodes.clear();
//         edges.clear();
//         first_parental_index = 0;
//     }
};

#endif
