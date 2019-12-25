#pragma once
// an original graph library
#include <string>

namespace simplex {

  template <typename iterator_t>
  class iter_pair {
    iterator_t _begin, _end;
  public:
    iter_pair(const iterator_t &begin, const iterator_t &end)
      : _begin{begin}, _end{end} {}
    const iterator_t &begin() const { return _begin; }
    const iterator_t &end() const { return _end; }
  };

  template <typename Derived>
  struct GraphTraits;

#define _G_THIS static_cast<Derived*>(this)
#define _G_CTHIS static_cast<const Derived*>(this)

  template <typename Derived>
  class _GraphCRTP {
  public:
    using vertex = typename GraphTraits<Derived>::vertex;
    using edge = typename GraphTraits<Derived>::edge;
    // all iters are const iters
    using vertex_iter = typename GraphTraits<Derived>::vertex_iter;
    using edge_iter = typename GraphTraits<Derived>::edge_iter;
    using adjacent_iter = typename GraphTraits<Derived>::adjacent_iter;
    using incident_iter = typename GraphTraits<Derived>::incident_iter;

    inline void
    clear()
    { _G_THIS->clear(); }

    inline iter_pair<vertex_iter>
    vertices() const
    { return _G_CTHIS->vertices(); };

    inline iter_pair<edge_iter>
    edges() const
    { return _G_CTHIS->edges(); };

    inline iter_pair<adjacent_iter>
    adjacent(const vertex &v) const
    { return _G_CTHIS->adjacent(v); }

    inline iter_pair<incident_iter>
    incident(const vertex &v) const
    { return _G_CTHIS->incident(v); }

    // add if does not exist, do nothing if it does
    void
    add_vertex(const vertex &v)
    { _G_THIS->add_vertex(v); }

    bool
    contains_vertex(const vertex &v) const
    { return _G_CTHIS->contains_vertex(v); }

    // will add the vertices if they do not exist
    void
    add_edge(const edge &e)
    { _G_THIS->add_edge(e); }

    bool
    contains_edge(const edge &e) const
    { return _G_CTHIS->contains_edge(e); }

    size_t
    vertex_to_id(const vertex &v) const
    { return _G_CTHIS->vertex_to_id(v); }

    const vertex&
    id_to_vertex(size_t id) const
    { return _G_CTHIS->id_to_vertex(id); }
  };

#undef _G_THIS
#undef _G_CTHIS
}


