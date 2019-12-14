#pragma once

#include "graph.hpp"
#include "../utils/disjoint_set.hpp"

#include <algorithm>

namespace simplex {

  // weighted kruskals
  template <typename Derived>
  _GraphCRTP<Derived>
  kruskals(const _GraphCRTP<Derived> &graph)
  {
    disjoint_set<typename _GraphCRTP<Derived>::vertex> ds;
    _GraphCRTP<Derived> mst;

    for (auto &v : graph.vertices())
      ds.make_set(v);
    auto list = graph.edges();
    std::sort(list.begin(), list.end(),
        [](auto &&e1, auto &&e2) { return e1.weight < e2.weight; });
    for (auto &e : list) {
      if (ds.find(e.first()) != ds.find(e.second())) {
        mst.add_edge(e);
        ds.union_set(e.first(), e.second());
      }
    }
    return mst;
  }

}

