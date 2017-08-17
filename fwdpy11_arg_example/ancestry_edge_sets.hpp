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

std::pair<std::vector<std::pair<double, double>>,
          std::vector<std::pair<double, double>>>
split_breakpoints(const std::vector<double>& breakpoints,
                  const double start = 0., const double stop = 1.);

bool
overlap(const std::tuple<double, double>& lhs,
        const std::tuple<double, double>& rhs);

#endif
