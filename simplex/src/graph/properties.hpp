#pragma once

#include "depth_first_search.hpp"

#include <unordered_set>

namespace simplex {
  template <typename Graph>
  size_t
  num_connected_components(const Graph &g) {
    struct _Visitor {
      std::unordered_set<typename Graph::vertex> seen;
      void discover_vertex(const typename Graph::vertex &v, const Graph &g) {
        seen.insert(v);
      }
    } visitor;
    size_t cnt = 0;
    for (auto &&v : g.vertices()) {
      if (visitor.seen.find(v) == visitor.seen.end()) {
        ++cnt;
        depth_first_search(g, v, visitor);
      }
    }
    return cnt;
  }
}

