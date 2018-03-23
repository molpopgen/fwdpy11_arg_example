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

inline void
reverse_time(std::vector<node>& nodes)
{
    if (nodes.empty())
        return;

    //convert forward time to backwards time
    auto max_gen = nodes.back().generation;

    for (auto& n : nodes)
        {
            n.generation -= max_gen;
            n.generation *= -1.0;
        }
}

struct ancestry_data
{
    using integer_type = decltype(edge::parent);
    /// Nodes:
    std::vector<node> nodes;
    /// The ARG:
    std::vector<edge> edges;
    /// Mutations:
    std::vector<mutation> mutations;
    std::vector<integer_type> samples;
    pybind11::object lock_;
    ancestry_data() : nodes{}, edges{}, mutations{}, samples{}, lock_{}
    {
        pybind11::module threading = pybind11::module::import("threading");
        lock_ = threading.attr("Lock")();
    }
};

struct ancestry_tracker
{
    using integer_type = decltype(edge::parent);
    /// Nodes:
    std::vector<node> nodes;
    /// The ARG:
    std::vector<edge> edges;
    /// The edges generated for each generation:
    std::vector<edge> temp;
    /// Mutations:
    std::vector<mutation> mutations;
    /// This is used as the sample indexes for msprime:
    std::vector<integer_type> offspring_indexes;
    integer_type generation, next_index, first_parental_index;
    std::uint32_t lastN;
    decltype(node::generation) last_gc_time;
    ancestry_tracker(const integer_type N, const integer_type next_index_)
        : nodes{ std::vector<node>() }, edges{ std::vector<edge>() },
          temp{ std::vector<edge>() }, mutations{ std::vector<mutation>() }, 
          offspring_indexes{ std::vector<integer_type>() }, generation{ 1 },
          next_index{ next_index_ }, first_parental_index{ 0 },
          lastN{ static_cast<std::uint32_t>(N) }, last_gc_time{ 0.0 }
    {
        nodes.reserve(2 * N);
        edges.reserve(2 * N);
        temp.reserve(N);
        //no need to reserve mutation space

        //Initialize 2N nodes for the generation 0
        //if next_index == 0, then did not initialize with tree sequence
        if (next_index_ == 0)
            {
                for (integer_type i = 0; i < 2 * N; ++i)
                    {
                        //ID, time 0, population 0
                        nodes.emplace_back(make_node(i, 0.0, 0));
                    }
            }
    }

    std::tuple<integer_type, integer_type>
    get_parent_ids(const std::uint32_t p, const int did_swap)
    {
        return std::make_tuple(
            first_parental_index + 2 * static_cast<integer_type>(p) + did_swap,
            first_parental_index + 2 * static_cast<integer_type>(p)
                + !did_swap);
    }

    std::tuple<integer_type, integer_type>
    get_next_indexes()
    {
        auto rv = std::make_tuple(next_index, next_index + 1);
        next_index += 2;
        offspring_indexes.push_back(std::get<0>(rv));
        offspring_indexes.push_back(std::get<1>(rv));
        return rv;
    }

    void
    add_edges(const std::vector<std::pair<double, double>>& breakpoints,
              const integer_type parent, const integer_type child)
    {
        for (auto&& bi : breakpoints)
            {
                temp.emplace_back(
                    make_edge(bi.first, bi.second, parent, child));
            }
    }

    void
    add_mutations(const std::vector<std::uint32_t>& mutation_ids, const integer_type node_id)
    {
        for (auto&& mut_id : mutation_ids)
            {
                mutations.emplace_back(
                    make_mutation(node_id, mut_id));
            }
    }

    void
    finish_generation()
    {
        for (auto&& oi : offspring_indexes)
            {
                nodes.emplace_back(make_node(oi, generation, 0));
            }
        edges.insert(edges.end(), temp.begin(), temp.end());
        lastN = next_index - first_parental_index;
        first_parental_index = offspring_indexes.front();

        temp.clear();
        ++generation;
    }

    void
    post_process_gc(pybind11::tuple t, const bool clear = true)
    {
        pybind11::bool_ gc = t[0].cast<bool>();
        if (!gc)
            return;

        last_gc_time = generation;
        next_index = t[1].cast<integer_type>();
        // establish last parental index:
        first_parental_index = 0;
        if (clear)
            {
                nodes.clear();
                edges.clear();
                mutations.clear();
            }
    }

    void
    exchange_for_async(ancestry_data& a)
    {
        nodes.swap(a.nodes);
        edges.swap(a.edges);
        mutations.swap(a.mutations);
        a.samples.assign(offspring_indexes.begin(), offspring_indexes.end());
        nodes.clear();
        edges.clear();
        first_parental_index = 0;
    }
};

#endif
