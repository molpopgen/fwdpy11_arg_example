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
#include <map>
#include <limits>
#include <cstdint>
#include <atomic>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

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
    using index_vec = std::vector<integer_type>;
    /// Nodes:
    std::vector<node> nodes;
    /// The ARG:
    std::vector<edge> edges;
    /// The edges generated for each generation:
    std::vector<edge> temp;
    /// Mutations:
    std::vector<mutation> mutations;
    /// strided per-population, per-generation indexes for argsimplifier:
    index_vec pop_gen_indexes;
    // current generation, total generation, next node ID to use, current index generation
    integer_type generation, total_generations, next_index, index_gen;
    ancestry_tracker(const integer_type N, 
                     const integer_type next_index_,
                     const integer_type total_generations_)
        : nodes{ std::vector<node>() }, edges{ std::vector<edge>() },
          temp{ std::vector<edge>() }, mutations{ std::vector<mutation>() },
          pop_gen_indexes{ index_vec(1) }, 
          generation{ 1 }, total_generations{ total_generations_ },
          next_index{ next_index_ }, index_gen{ 1 }
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
        pop_gen_indexes.emplace_back(next_index);
    }

    std::tuple<integer_type, integer_type>
    get_parent_ids(const std::uint32_t p, const int did_swap)
    {
        integer_type first_parental_index = pop_gen_indexes[index_gen-1];
        return std::make_tuple(
            first_parental_index + 2 * static_cast<integer_type>(p) + did_swap,
            first_parental_index + 2 * static_cast<integer_type>(p)
                + !did_swap);
    }

    std::tuple<integer_type, integer_type>
    get_next_indexes()
    {
        auto rv = std::make_tuple(next_index, next_index + 1);
        double rev_gen = total_generations - generation;
        nodes.emplace_back(node{std::get<0>(rv), 0, rev_gen});
        nodes.emplace_back(node{std::get<1>(rv), 0, rev_gen});
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

    void
    add_mutations(const std::vector<std::uint32_t>& mutation_ids, const integer_type node_id)
    {
        for (auto&& mut_id : mutation_ids)
            {
                mutations.emplace_back(mutation{node_id, mut_id});
            }
    }

    void
    finish_generation()
    {
        edges.insert(edges.end(), temp.begin(), temp.end());
        temp.clear();
        pop_gen_indexes.emplace_back(next_index);
        ++index_gen;
        ++generation;
    }

    void
    post_process_gc(pybind11::tuple t, const bool clear = true)
    {
        pybind11::bool_ gc = t[0].cast<bool>();
        if (!gc)
            return;

        next_index = t[1].cast<integer_type>();
        // reset last parental index to 0 (doesn't matter for earlier indices):
        // this is the last generation to be processed by GC 
        index_gen = 1;
        // reset current parental index to next_index:
        pop_gen_indexes.resize(1);
        pop_gen_indexes.emplace_back(next_index);
        
        if (clear)
            {
                nodes.clear();
                edges.clear();
                mutations.clear();
            }
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
