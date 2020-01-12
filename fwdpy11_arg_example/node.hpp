#ifndef ANCESTRY_NODE_HPP__
#define ANCESTRY_NODE_HPP__

#include <cstdint>
#include <tuple>

struct node
{
    std::int32_t id;
	std::int32_t population;
	double generation;
};

inline auto
get_tied_node(const node& n) -> decltype(std::tie(n.generation, n.population, n.id))
{
    return std::tie(n.generation, n.id, n.population);
}

inline bool
operator<(const node& lhs, const node& rhs)
{
    // sort order is generation, then population, then id
    return get_tied_node(lhs) < get_tied_node(rhs);
}

#endif
