#ifndef ANCESTRY_NODE_HPP__
#define ANCESTRY_NODE_HPP__

#include <cstdint>
#include <tuple>

struct node
{
    std::uint32_t id, generation;
    std::uint32_t deme;
};

inline node
make_node(std::uint32_t id, std::uint32_t generation, std::uint32_t deme)
{
    node n;
    n.id = id;
    n.generation = generation;
    n.deme = deme;
    return n;
}

inline auto
get_tied_node(const node& n) -> decltype(std::tie(n.generation, n.deme, n.id))
{
    return std::tie(n.generation, n.id, n.deme);
}

inline bool
operator<(const node& lhs, const node& rhs)
{
    // sort order is generation, then deme, then id
    return get_tied_node(lhs) < get_tied_node(rhs);
}

#endif
