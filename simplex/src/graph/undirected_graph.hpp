#pragma once

#include <iostream>
#include "base_graph.hpp"
#include "adjacency_list.hpp"

namespace simplex {

  // turns an directed graph impl to undirected
  // the directed graph edge iterator must have an order
  template <typename GraphImpl>
  class Undirected: public _GraphCRTP<GraphImpl>
  {
    using _this_t = Undirected<GraphImpl>;
    GraphImpl _impl;
  public:
    using vertex = typename GraphTraits<_this_t>::vertex;
    using edge = typename GraphTraits<_this_t>::edge;
    using vertex_iter = typename GraphTraits<_this_t>::vertex_iter;
    using edge_iter = typename GraphTraits<_this_t>::edge_iter;
    using adjacent_iter = typename GraphTraits<_this_t>::adjacent_iter;
    using incident_iter = typename GraphTraits<_this_t>::incident_iter;

    inline void clear() { _impl.clear(); }

    inline iter_pair<vertex_iter>
    vertices() const 
    { return _impl.vertices(); }

    iter_pair<edge_iter>
    edges() const
    {
      auto _e = _impl.edges();
      return iter_pair<edge_iter>{edge_iter{_e.begin(), _impl},
                                  edge_iter{_e.end(), _impl}};
    }

    inline iter_pair<adjacent_iter>
    adjacent(const vertex &v) const
    { return _impl.adjacent(v); }

    inline iter_pair<incident_iter>
    incident(const vertex &v) const
    { return _impl.incident(v); }

    inline void
    add_vertex(const vertex &v)
    { _impl.add_vertex(v); }

    inline bool
    contains_vertex(const vertex &v)
    { return _impl.contains_vertex(v); }

    void
    add_edge(const edge &e)
    {
      if (contains_edge(e)) return;
      _impl.add_edge(e);
      _impl.add_edge(edge{e.second, e.first, e.weight});
    }

    inline bool
    contains_edge(const edge &e)
    { return _impl.contains_edge(e); }

    size_t
    vertex_to_id(const vertex &v) const
    { return _impl.vertex_to_id(v); }

    const vertex&
    id_to_vertex(size_t id) const
    { return _impl.id_to_vertex(id); }

    // TODO use template metaprogramming to disable these default implementations
    // and use the derived implementations, if they exist
    size_t
    order() const
    {
      size_t cnt = 0;
      auto vs = vertices();
      for (auto it = vs.begin(); it != vs.end(); ++it) ++cnt;
      return cnt;
    }

    inline bool
    empty() const
    {
      auto v = vertices();
      if (v.begin() == v.end()) return true;
      return false;
    }

    size_t
    size() const
    {
      size_t cnt = 0;
      for (auto &&e : edges()) ++cnt;
      return cnt;
    }

    friend std::ostream&
    operator<<(std::ostream &o, const _this_t &graph)
    { return operator<<(o, graph.impl); }
  };

  template <typename GraphImpl>
  struct GraphTraits<Undirected<GraphImpl>> {
  private:
    class _edge_iter;
  public:
    using vertex = typename GraphImpl::vertex;
    using edge = typename GraphImpl::edge;
    using vertex_iter = typename GraphImpl::vertex_iter;
    using edge_iter = _edge_iter;
    using adjacent_iter = typename GraphImpl::adjacent_iter;
    using incident_iter = typename GraphImpl::incident_iter;
  private:
    class _edge_iter{
      using _edge_iter_t = typename GraphImpl::edge_iter;
      _edge_iter_t it;
      const GraphImpl *_graph;
    public:
      explicit _edge_iter(const _edge_iter_t &it, const GraphImpl &graph)
        : it{it}, _graph{&graph}
      {}

      _edge_iter(const _edge_iter &other)
        : it{other.it}, _graph{other._graph}
      {}

      bool
      operator!=(const _edge_iter &other) const
      { return it != other.it; }

      _edge_iter&
      operator++()
      {
        auto e = _graph->edges().end();
        do {
          ++it;
          if (!(it != e)) return *this;
        } while (_graph->vertex_to_id((*it).first)
               > _graph->vertex_to_id((*it).second));
        return *this;
      }

      edge
      operator*()
      { return *it; }
    };
  };
}

