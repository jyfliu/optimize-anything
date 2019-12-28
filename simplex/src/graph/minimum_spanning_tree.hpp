#pragma once

#include "base_graph.hpp"
#include "../utils/disjoint_set.hpp"

#include <vector>
#include <algorithm>

namespace simplex {

  // weighted kruskals
  template <typename Derived>
  Derived
  kruskals(const Derived &graph)
  {
    disjoint_set<typename _GraphCRTP<Derived>::vertex> ds;
    Derived mst;

    for (auto &&v : graph.vertices()) {
      ds.make_set(v);
    }
    auto itp = graph.edges();
    // TODO add a graph.num_edges and reserve that space in the vector
    auto list = std::vector<typename Derived::edge>{};
    for (auto &&e : itp) {
      list.emplace_back(e);
    }
    std::sort(list.begin(), list.end(),
        [](auto &&e1, auto &&e2) { return e1.weight < e2.weight; });
    for (auto &&e : list) {
      if (ds.find(e.first) != ds.find(e.second)) {
        mst.add_edge(e);
        ds.union_set(e.first, e.second);
      }
    }
    return mst;
  }

}

