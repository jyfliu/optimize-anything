#pragma once

#include "base_graph.hpp"
#include "depth_first_search.hpp"

#include <unordered_set>

namespace simplex {

  class not_bipartite{};

  template <typename GraphImpl>
  class Bipartition: public _GraphCRTP<Bipartition<GraphImpl>> {
    GraphImpl _impl;
    std::unordered_set<typename GraphImpl::vertex> _w, _j;
    struct _Visitor {
      std::unordered_set<typename GraphImpl::vertex> w, j;
      bool fail = false;
      void tree_edge(const typename GraphImpl::edge &e, const Bipartition &g) {
        if (w.find(e.first) != w.end()) {
          if (w.find(e.second) != w.end()) throw not_bipartite{};
          j.insert(e.second);
        } else if (j.find(e.first) != j.end()) {
          if (j.find(e.second) != j.end()) throw not_bipartite{};
          w.insert(e.second);
        }
        else throw; // shouldn't happen
      }
    };
  public:
    using vertex = typename GraphImpl::vertex;
    using edge = typename GraphImpl::edge;
    using vertex_iter = typename GraphImpl::vertex_iter;
    using edge_iter = typename GraphImpl::edge_iter;
    using adjacent_iter = typename GraphImpl::adjacent_iter;
    using incident_iter = typename GraphImpl::incident_iter;

    explicit
    Bipartition(const GraphImpl &derived)
      : _impl{derived}
    { _colour(); }

    inline void
    clear()
    { _impl.clear(); }

    inline iter_pair<vertex_iter>
    vertices() const
    { return _impl.vertices(); }

    inline iter_pair<edge_iter>
    edges() const
    { return _impl.edges(); }

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
    contains_vertex(const vertex &v) const
    { return _impl.contains_vertex(v); }

    inline void
    add_edge(const edge &e)
    { _impl.add_edge(e); }

    inline bool
    contains_edge(const edge &e) const
    { return _impl.contains_edge(e); }

    inline size_t
    vertex_to_id(const vertex &v) const
    { return _impl.vertex_to_id(v); }

    inline const vertex&
    id_to_vertex(size_t id) const
    { return _impl.id_to_vertex(id); }

    inline size_t
    order() const
    { return _impl.order(); }

    inline bool
    empty() const
    { return _impl.empty(); }

    inline size_t
    size() const
    { return _impl.size(); }

    const std::unordered_set<vertex> &
    W() const
    { return _w; }

    const std::unordered_set<vertex> &
    J() const
    { return _j; }
  private:
    void
    _colour()
    {
      _Visitor visitor;
      for (auto &&v : _impl.vertices()) {
        if (visitor.w.find(v) == visitor.w.end()
            and visitor.j.find(v) == visitor.j.end()) {
          if (visitor.w.size() < visitor.j.size()) {
            visitor.w.insert(v);
          } else {
            visitor.j.insert(v);
          }
          depth_first_search(*this, v, visitor);
        }
      }
      _w = std::move(visitor.w);
      _j = std::move(visitor.j);
      return;
    }
  };

  template <typename GraphImpl>
  struct GraphTraits<Bipartition<GraphImpl>> {
    using vertex = typename GraphImpl::vertex;
    using edge = typename GraphImpl::edge;
    using vertex_iter = typename GraphImpl::vertex_iter;
    using edge_iter = typename GraphImpl::edge_iter;
    using adjacent_iter = typename GraphImpl::adjacent_iter;
    using incident_iter = typename GraphImpl::incident_iter;
  };

}
