#ifndef ANCESTRY_ANCESTRY_TRACKER_HPP__
#define ANCESTRY_ANCESTRY_TRACKER_HPP__

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <pybind11/pybind11.h>

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
    integer_type generation, next_index, first_parental_index,
        first_child_index;
    std::uint32_t lastN;
    ancestry_tracker(const integer_type N)
        : nodes{ std::vector<node>() }, edges{ std::vector<edge>() },
          temp{ std::vector<edge>() },
          //parental_indexes{ std::vector<integer_type>() },
          offspring_indexes{ std::vector<integer_type>() }, generation{ 1 },
          lastN{ N }, next_index{ 2 * N }, first_parental_index{ 0 },
          first_child_index{ 2 * N }
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
        // for(auto && i : reproduced)
        //{
        //	std::cout << i << ' ';
        //}
        // std::cout << " -> ";
        // std::cerr << "G = " << generation << '\n';
        // prange
        //    = std::equal_range(edges.begin(), edges.end(), temp.front().pgen,
        //                       get_parental_generation());
        // edges.erase(std::remove_if(prange.first, prange.second,
        //                           [this](const edge& e) {
        //                               return reproduced.find(e.cid)
        //                                      == reproduced.end();
        //                           }),
        //            prange.second);
        // std::cout << generation << ' '<< edges.size() << ' ' <<
        // std::distance(prange.first,prange.second) << '\n';
        // integer_type extinct=0;
        // std::cout << reproduced.size() << " reproduced, " <<
        //	std::distance(prange.first,prange.second) << '\n';
        // for(auto i = prange.first;i!=prange.second;++i)
        //{
        //	if(reproduced.find(i->cid) == reproduced.end())
        //	{
        //		++extinct;
        //	}
        //}
        // std::cout << extinct << " extinct lineages\n";
        //add_nodes();

        for (auto&& oi : offspring_indexes)
            {
                nodes.emplace_back(make_node(oi, generation, 0));
            }
        std::sort(temp.begin(), temp.end());
        edges.insert(edges.end(), temp.begin(), temp.end());
        lastN = next_index - first_parental_index;
        first_parental_index = first_child_index;
        first_child_index = next_index;

        //parental_indexes.swap(offspring_indexes);
        //offspring_indexes.clear();
        temp.clear();
        ++generation;
    }

    void
    prep_for_gc()
    {
		//Sorting the nodes is easy.
		//For the edges, I should sort on parent id in
		//decreasing order and then increasing child id
		//within that.  That means there's not a "natural"
		//operator> or < that would work, and I should farm
		//it off to a new struct.
        if (nodes.empty())
            return;

        auto max_gen = nodes.back().generation;
        //std::sort(nodes.begin(), nodes.end());
        for (auto& n : nodes)
            {
                n.generation -= max_gen;
                n.generation *= -1.0;
            }
        //std::sort(edges.begin(), edges.end(),
        //          [](const edge& lhs, const edge& rhs) { return lhs > rhs; });
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
        std::sort(temp.begin(), temp.end());
        std::vector<std::tuple<double, double>> rv;
        for (auto i = temp.begin() + 1; i != temp.end(); ++i)
            {
                rv.emplace_back(*(i - 1), *i);
            }
        return rv;
    }
};

#endif
