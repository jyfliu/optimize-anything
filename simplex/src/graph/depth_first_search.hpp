#pragma once

#include "base_graph.hpp"

#include <vector>
#include <stack>

#ifndef dfs_USING_IF_EXIST
// use expression sfinae to conditionally include functions
#define dfs_USING_IF_EXIST(fn_name, param_t)\
  private:\
    template <typename X=Visitor>\
    inline auto\
    _ ## fn_name(const param_t &_p, const Graph &_g, int)\
      -> decltype(static_cast<X&>(_visitor).fn_name(_p, _g), void())\
    { _visitor.fn_name(_p, _g); }\
    \
    template <typename X=Visitor>\
    inline void\
    _ ## fn_name(const param_t &, const Graph &, long)\
    {}\
    \
  public:\
    inline void\
    fn_name(const param_t &_p, const Graph &_g)\
    { _ ## fn_name(_p, _g, 0); }

namespace simplex {
  // wraps a visitor to provide potentially unimplemented functions
  template <typename Visitor, typename Graph>
  class _DFSVisitor {
    using vertex = typename Graph::vertex;
    using edge = typename Graph::edge;
    Visitor &_visitor;
  public:
    explicit _DFSVisitor(Visitor &visitor): _visitor{visitor} {}
    _DFSVisitor &operator=(_DFSVisitor) = delete;

    // boost options
    dfs_USING_IF_EXIST(init_vertex, vertex);
    dfs_USING_IF_EXIST(start_vertex, vertex);
    dfs_USING_IF_EXIST(discover_vertex, vertex);
    dfs_USING_IF_EXIST(examine_edge, edge);
    dfs_USING_IF_EXIST(tree_edge, edge);
    dfs_USING_IF_EXIST(back_edge, edge);
    dfs_USING_IF_EXIST(forward_or_cross_edge, edge);
    dfs_USING_IF_EXIST(finish_vertex, vertex);

    Visitor &expose() const { return _visitor; }
  };

  enum _DFSColors { WHITE, GRAY, BLACK };

  template <typename Derived,
            typename Visitor,
            typename Stack=std::stack<typename Derived::vertex,
                                      std::vector<typename Derived::vertex>>>
  Visitor
  depth_first_search(const Derived &graph,
                     const typename Derived::vertex &start,
                     Visitor &_visitor)
  {
    using vertex = typename Derived::vertex;
    size_t order = graph.order();
    Stack stack;
    std::vector<_DFSColors> visited(order);
    _DFSVisitor<Visitor, Derived> visitor(_visitor);

    stack.push(start);
    while (!stack.empty()) {
      vertex u = stack.top();
      stack.pop();
      visited[graph.vertex_to_id(u)] = GRAY;                                    visitor.discover_vertex(u, graph);
      for (auto &&e : graph.incident(u)) {
        vertex v = e.second;                                                    visitor.examine_edge(e, graph);
        switch(visited[graph.vertex_to_id(v)]) {
          case WHITE:                                                           visitor.tree_edge(e, graph);
            stack.push(v);
            break;
          case GRAY:                                                            visitor.back_edge(e, graph);
            break;
          case BLACK:                                                           visitor.forward_or_cross_edge(e, graph);
            break;
        }
      }
      visited[graph.vertex_to_id(u)] = BLACK;                                   visitor.finish_vertex(u, graph);
    }
    return visitor.expose();
  }
}
#endif
#undef dfs_USING_IF_EXIST

