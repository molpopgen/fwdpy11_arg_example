#ifndef ANCESTRY_EDGE_SETS_HPP__
#define ANCESTRY_EDGE_SETS_HPP__

#include <cstdint>
#include <tuple>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>
#include <cassert>

#include "node.hpp"
#include "edge.hpp"
#include "ancestry_tracker.hpp"

inline
std::pair<std::vector<std::pair<double, double>>,
          std::vector<std::pair<double, double>>>
split_breakpoints(const std::vector<double>& breakpoints,
                  const double start = 0., const double stop = 1.)
{
    std::vector<std::pair<double, double>> r1(
        1, std::make_pair(start, breakpoints[0])),
        r2;
    for (unsigned j = 1; j < breakpoints.size(); ++j)
        {
            double a = breakpoints[j - 1];
            double b = (j < breakpoints.size() - 1) ? breakpoints[j] : stop;
            if (j % 2 == 0.)
                {
                    r1.emplace_back(a, b);
                }
            else
                {
                    r2.emplace_back(a, b);
                }
        }
    return std::make_pair(std::move(r1), std::move(r2));
}

inline bool
overlap_detail(const std::tuple<double, double>& lhs,
               const std::tuple<double, double>& rhs)
{
    // return not (x[1] < y[0] or y[1] < x[0])
    return !((std::get<1>(lhs) <= std::get<0>(rhs))
             || (std::get<1>(rhs) <= std::get<0>(lhs)));
}

inline bool
overlap(const std::tuple<double, double>& lhs,
        const std::tuple<double, double>& rhs)
{
    if (std::get<0>(lhs) < std::get<0>(rhs))
        {
            return overlap_detail(lhs, rhs);
        }
    return overlap_detail(rhs, lhs);
}


inline
std::vector<std::tuple<double, double>>
sorted_tree_edges(const std::vector<edge>& edges)
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

#endif
