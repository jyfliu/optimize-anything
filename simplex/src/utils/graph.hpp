#pragma once
// an original graph library
#include <string>

namespace simplex {
  template <typename Label_t=int>
  struct Vertex {
    Label_t label;
  };

  template <typename Vertex_t=Vertex<>>
  struct VertexIterator {
    bool operator!=(const VertexIterator &other);
    VertexIterator &operator++();
    Vertex_t &operator*();
  };

  template <typename VertexIterator_t=VertexIterator<>>
  struct IteratorPair {
    const VertexIterator_t &begin() const;
    const VertexIterator_t &end() const;
  };

  template <typename Vertex_t=Vertex<>>
  struct Edge {
    Vertex_t first, second;
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

    template <typename iterator_t>
    struct iter_pair {
      const iterator_t &begin() const;
      const iterator_t &end() const;
    };

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
    { _G_THIS->contains_vertex(v); }

    // will add the vertices if they do not exist
    void
    add_edge(const edge &e)
    { _G_THIS->add_edge(e); }
  };

#undef _G_THIS
}

