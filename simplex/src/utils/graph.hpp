#pragma once
// an original graph library
#include <string>

namespace simplex {

  template <typename iterator_t>
  class iter_pair {
    const iterator_t &_begin, &_end;
  public:
    iter_pair(const iterator_t &begin, const iterator_t &end)
      : _begin{begin}, _end{end} {}
    const iterator_t &begin() const;
    const iterator_t &end() const;
  };

#define _G_THIS static_cast<Derived*>(this)

  template <typename Derived>
  class _GraphCRTP {
  public:
    using vertex = typename Derived::vertex;
    using edge = typename Derived::edge;
    // all iters are const iters
    using vertex_iter = typename Derived::vertex_iter;
    using edge_iter = typename Derived::edge_iter;
    using adjacent_iter = typename Derived::adjacent_iter;
    using incident_iter = typename Derived::incident_iter;

    inline void
    clear()
    { _G_THIS->clear(); }

    inline iter_pair<vertex_iter>
    vertices() const
    { return _G_THIS->vertices(); };

    inline iter_pair<edge_iter>
    edges() const
    { return _G_THIS->edges(); };

    inline iter_pair<adjacent_iter>
    adjacent(const vertex &v) const
    { return _G_THIS->adjacent(v); }

    inline iter_pair<incident_iter>
    incident(const vertex &v) const
    { return _G_THIS->incident(v); }

    // add if does not exist, do nothing if it does
    void
    add_vertex(const vertex &v)
    { _G_THIS->add_vertex(v); }

    bool
    contains_vertex(const vertex &v)
    { return _G_THIS->contains_vertex(v); }

    // will add the vertices if they do not exist
    void
    add_edge(const edge &e)
    { _G_THIS->add_edge(e); }

    bool
    contains_edge(const edge &e)
    { return _G_THIS->contains_edge(e); }
  };

#undef _G_THIS
}


