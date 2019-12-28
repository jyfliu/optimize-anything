#pragma once

#include "base_graph.hpp"

#include <unordered_set>

namespace simplex {

  template <typename Derived>
  class Bipartition: public Derived {
    std::unordered_set<typename Derived::vertex> _w, _j;
    using vertex = typename Derived::vertex;
  public:
    explicit
    Bipartition(const Derived &derived)
    {
      for (auto &&v : derived.vertices())
      { add_vertex(v); }
      for (auto &&e : derived.edges())
      { add_edge(e); }
      _colour();
    }
  private:
    void
    _colour()
    {
      auto v = Derived::vertices();
      if (v.begin() == v.end()) return;
      auto start = *v.begin();
    }
  };

}
