#ifndef ANCESTRY_ANCESTRY_TRACKER_HPP__
#define ANCESTRY_ANCESTRY_TRACKER_HPP__

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <cstdint>

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
    std::vector<integer_type> parental_indexes, offspring_indexes;
    std::pair<std::vector<edge>::iterator, std::vector<edge>::iterator> prange;
    integer_type generation, lastN, next_index, last_index_used;
    ancestry_tracker(const integer_type N)
        : nodes{ std::vector<node>() }, edges{ std::vector<edge>() },
          temp{ std::vector<edge>() },
          parental_indexes{ std::vector<integer_type>() },
          offspring_indexes{ std::vector<integer_type>() }, generation{ 1 },
          lastN{ N }, next_index{ 2 * N }, last_index_used{ 0 }
    {
        nodes.reserve(2 * N);
        edges.reserve(2 * N);
        temp.reserve(N);
		
		//Initialize 2N nodes for the generation 0
		for(integer_type i = 0 ; i < 2*N ; ++i)
		{
			//ID, time 0, population 0
			nodes.emplace_back(make_node(i, 0.0, 0));
		}
    }

    void
    initialize_generation()
    {
        last_index_used = next_index;
    }

    std::tuple<integer_type, integer_type>
    get_parent_ids(const std::uint32_t p, const int did_swap)
    {
        if (parental_indexes.empty())
            {
                return std::make_tuple(2 * p + did_swap, 2 * p + !did_swap);
            }
        return std::make_tuple(parental_indexes[2 * p + did_swap],
                               parental_indexes[2 * p + !did_swap]);
    }

    std::tuple<integer_type, integer_type>
    get_next_indexes(const bool update_offspring_indexes = true)
    {
        auto rv = std::make_tuple(next_index, next_index + 1);
        next_index += 2;
        if (update_offspring_indexes)
            {
                offspring_indexes.push_back(std::get<0>(rv));
                offspring_indexes.push_back(std::get<1>(rv));
            }
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
		for(auto && oi : offspring_indexes)
		{
			nodes.emplace_back(make_node(oi,generation,0));
		}
        std::sort(temp.begin(), temp.end());
        edges.insert(edges.end(), temp.begin(), temp.end());
        parental_indexes.swap(offspring_indexes);
        offspring_indexes.clear();
        temp.clear();
        ++generation;
    }

    void
    reconcile_for_msprime()
    {
        std::cout << "some checking: " << generation << ' ' << last_index_used
                  << '\n';
        std::sort(nodes.begin(), nodes.end());
        auto last_size = nodes.size();
        auto last_gen = generation - 1;
        std::unordered_set<integer_type> used;
        for (auto i = edges.rbegin(); i != edges.rend(); ++i)
            {
                if (i->child >= last_index_used)
                    {
                        if (used.find(i->child) == used.end())
                            {
                                nodes.emplace_back(
                                    make_node(i->child, last_gen, 0));
                                used.insert(i->child);
                            }
                    }
                else
                    {
                        break;
                    }
            }
        std::sort(nodes.begin() + last_size, nodes.end());
        std::cout << last_gen << '\n';
        // remove all edges where child IDs are not listed as parents in nodes
        // and the node time is < last_gen
        std::cout << "Number of edges: " << edges.size() << " -> ";
        edges.erase(std::remove_if(
                        edges.begin(), edges.end(),
                        [this, last_gen](const edge& e) {
                            //std::cout << e.parent << ' ' << e.child << ": ";
                            auto x = std::lower_bound(
                                nodes.begin(), nodes.end(), e.child,
                                [](const node& n, const integer_type c) {
                                    return n.id < c;
                                });
                            //if(x != nodes.end())
                            //{
                            //std::cout << x->id << ' ' << x->generation << ' ' << last_gen << ' ';
                            //}
                            if (x == nodes.end() || (x->generation < last_gen
                                                     && x->id != e.child))
                                {
                                    //        std::cout << "true\n";
                                    return true;
                                }
                            //std::cout << "false\n";
                            return false;
                        }),
                    edges.end());
        std::unordered_set<integer_type> used_ids_in_edges;
        for (auto&& e : edges)
            {
                used_ids_in_edges.insert(e.child);
                used_ids_in_edges.insert(e.parent);
            }
        std::cout << edges.size() << ", and " << nodes.size() << " nodes, and "
                  << used_ids_in_edges.size() << " unique ids in edges\n";
        std::cout << "these nodes are not in edges:\n";
        nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
                                   [&used_ids_in_edges](const node& n) {
                                       return used_ids_in_edges.find(n.id)
                                              == used_ids_in_edges.end();
                                   }),
                    nodes.end());
        std::cout << "Down to " << nodes.size() << " nodes\n";
        //for(auto && n : nodes)
        //{
        //	if(used_ids_in_edges.find(n.id) == used_ids_in_edges.end())
        //	{
        //		std::cout << n.id << ' ' << n.generation << '\n';
        //	}
        //}
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
