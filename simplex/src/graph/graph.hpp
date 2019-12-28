#pragma once

#include "base_graph.hpp"
#include "adjacency_list.hpp"
#include "undirected_graph.hpp"
#include "bipartite.hpp"

namespace simplex {
  template <typename Vertex_t>
  using Graph = Undirected<AdjacencyList<Vertex_t>>;
}

