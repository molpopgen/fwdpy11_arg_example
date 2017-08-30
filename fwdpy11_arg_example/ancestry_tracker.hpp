#ifndef ANCESTRY_ANCESTRY_TRACKER_HPP__
#define ANCESTRY_ANCESTRY_TRACKER_HPP__

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <limits>
#include <unordered_set>
#include <cstdint>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "node.hpp"
#include "edge.hpp"

struct ancestry_tracker
{
    using integer_type = decltype(edge::parent);
    std::vector<node> nodes;
    /// The ARG
    std::vector<edge> edges;
    /// Just in case.
    std::vector<edge> temp;
    //std::vector<integer_type> parental_indexes, offspring_indexes;
    std::vector<integer_type> offspring_indexes;
    std::pair<std::vector<edge>::iterator, std::vector<edge>::iterator> prange;
    integer_type generation, next_index, first_parental_index;
    std::uint32_t lastN;
    decltype(node::generation) last_gc_time;
    //std::unordered_map<integer_type, integer_type> sample_map;
    ancestry_tracker(const integer_type N)
        : nodes{ std::vector<node>() }, edges{ std::vector<edge>() },
          temp{ std::vector<edge>() },
          //parental_indexes{ std::vector<integer_type>() },
          offspring_indexes{ std::vector<integer_type>() }, generation{ 1 },
          next_index{ 2 * N }, first_parental_index{ 0 },
 lastN{ static_cast<std::uint32_t>(N) },
          last_gc_time{ 0.0 }
    {
        nodes.reserve(2 * N);
        edges.reserve(2 * N);
        temp.reserve(N);

        //Initialize 2N nodes for the generation 0
        for (integer_type i = 0; i < 2 * N; ++i)
            {
                //ID, time 0, population 0
                nodes.emplace_back(make_node(i, 0.0, 0));
            }
    }

    void
    initialize_generation()
    {
        //last_index_used = next_index;
    }

    std::tuple<integer_type, integer_type>
    get_parent_ids(const std::uint32_t p, const int did_swap)
    {
        return std::make_tuple(
            first_parental_index + 2 * static_cast<integer_type>(p) + did_swap,
            first_parental_index + 2 * static_cast<integer_type>(p)
                + !did_swap);
        // if (parental_indexes.empty())
        //     {
        //         return std::make_tuple(2 * p + did_swap, 2 * p + !did_swap);
        //     }
        // return std::make_tuple(parental_indexes[2 * p + did_swap],
        //                        parental_indexes[2 * p + !did_swap]);
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

    // void
    // add_nodes()
    // {
    //     std::unordered_set<integer_type> used;
    //     for (auto& e : temp)
    //         {
    //             if (used.find(e.parent) == used.end())
    //                 {
    //                     nodes.emplace_back(
    //                         make_node(e.parent, generation - 1, 0));
    //                     used.insert(e.parent);
    //                 }
    //         }
    // }

    void
    finish_generation()
    {
        for (auto&& oi : offspring_indexes)
            {
                nodes.emplace_back(make_node(oi, generation, 0));
            }
        //std::sort(temp.begin(), temp.end());
        edges.insert(edges.end(), temp.begin(), temp.end());
        lastN = next_index - first_parental_index;
		//pybind11::print("changing indexes from: ", first_parental_index, next_index);
        first_parental_index = offspring_indexes.front(); 
		//pybind11::print("changing indexes to: ", first_parental_index, next_index);

        //parental_indexes.swap(offspring_indexes);
        //offspring_indexes.clear();
        temp.clear();
        ++generation;
    }

    std::vector<node>
    update_nodes(pybind11::array_t<double> generations_from_msprime)
    //We need to get node info back from msprime
    //and update some data, like the node times.
    //more generally, we need to copy population, too,
    //but this is not a general updater...
    {
        auto g = generations_from_msprime.unchecked<1>();
        std::vector<node> tnodes;
        tnodes.reserve(g.shape(0));
        double delta = static_cast<double>(generation) - last_gc_time;
        for (std::size_t i = 0; i < g.shape(0); ++i)
            {
                //pybind11::print(g(i),last_gc_time,delta);
                tnodes.emplace_back(make_node(i, g(i) + delta, 0));
            }
        return tnodes;
    }

    void
    prep_for_gc()
    {
        //Sorting the nodes is easy.
        //To sort edges, we need to add parental
        //time to edge struct and update functions
        //accordingly.  Then, we can define a closure
        //that will give same results as msprime.sort_tables().
        //We may go this route to improve performance.
        if (nodes.empty())
            return;

        auto max_gen = nodes.back().generation;
        //std::sort(nodes.begin(), nodes.end());
        for (auto& n : nodes)
            {
                n.generation -= max_gen;
                n.generation *= -1.0;
            }
        //auto tnodes = update_nodes(generations_from_msprime);
        //auto x = tnodes.size();
        //tnodes.insert(tnodes.end(), std::make_move_iterator(nodes.begin()),
        //              std::make_move_iterator(nodes.end()));
        //if (x)
        //    {
        //        for (auto&& n : tnodes)
        //            {
        //                pybind11::print(n.id, n.generation);
        //            }
        //    }
        //nodes.swap(tnodes);
        //tnodes.clear();
        //std::sort(edges.begin(), edges.end(),
        //          [](const edge& lhs, const edge& rhs) {
        //              if (lhs.parent != rhs.parent)
        //                  {
        //                      return lhs.parent > rhs.parent;
        //                  }
        //              if (lhs.child != rhs.child)
        //                  {
        //                      return lhs.child > rhs.child;
        //                  }
        //			  return true;
        //              //return lhs.right < rhs.right;
        //          });
    }

    void
    post_process_gc(pybind11::tuple t)
    {
        pybind11::bool_ gc = t[0].cast<bool>();
        if (!gc)
            return;

        last_gc_time = generation;
		// pybind11::print("post_process_gc: ", generation);
        next_index = t[1].cast<integer_type>();
        //sample_map = t[2].cast<decltype(sample_map)>();
        // establish last parental index:
		first_parental_index = 0;
        nodes.clear();
        edges.clear();
    }

    std::vector<std::tuple<double, double>>
    sorted_tree_edges(const std::vector<edge>& edges)
    /// Not sure where this should be, so it is here for now
    {
        // Note: not general.
        // Should have a more flexible left boundary
        // than this
        std::unordered_set<double> intervals{ 0.0 };
        for (auto&& e : edges)
            {
                intervals.emplace(e.right);
            }
        std::vector<double> temp(intervals.begin(), intervals.end());
        //std::sort(temp.begin(), temp.end());
        std::vector<std::tuple<double, double>> rv;
        for (auto i = temp.begin() + 1; i != temp.end(); ++i)
            {
                rv.emplace_back(*(i - 1), *i);
            }
        return rv;
    }
};

#endif
