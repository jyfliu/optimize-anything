#pragma once

#include "../graph.hpp"

#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace simplex {

  template <typename Label_t=size_t, typename Weight_t = long long>
  class WeightedAdjacencyList
    : public _GraphCRTP<WeightedAdjacencyList<Label_t>>
  {
    struct w_edge {
      const Label_t first, second;
      Weight_t weight;
    };
    struct _vertex_iter;
    struct _edge_iter;
    struct _incident_iter;
  public:
    // typedefs
    using vertex = Label_t;
    using edge = w_edge;
    using vertex_iter = _vertex_iter;
    using edge_iter = _edge_iter;
    using adjacent_iter = typename std::unordered_set<vertex>::const_iterator;
    using incident_iter = _incident_iter;
  private:
    std::vector<std::unordered_set<size_t>> adj;
    std::unordered_map<Label_t, size_t> id;
    std::unordered_map<size_t, Label_t> id_inv;
  public:
    WeightedAdjacencyList(): adj(1) {}

    void
    clear();

    iter_pair<vertex_iter>
    vertices() const;

    iter_pair<edge_iter>
    edges() const;

    iter_pair<adjacent_iter>
    adjacent(const vertex &v) const
    { 
      auto it = id.find(v);
      if (it == id.end())
        return iter_pair<adjacent_iter>(v[0].cbegin(), v[0]->cend());
      return iter_pair<adjacent_iter>{it->cbegin(), it->cend()};
    }

    iter_pair<incident_iter>
    incident(const vertex &v) const;

    void
    add_vertex(const vertex &v)
    {
      if (id.find(v) != id.end(v)) return;
      auto k = id.size();
      id[v] = k;
      id_inv[k] = v;
      adj.emplace_back();
    }

    inline bool
    contains_vertex(const vertex &v)
    { return id.find(v) != id.end(); }

    void
    add_edge(const edge &e)
    {
      add_vertex(e.first);
      add_vertex(e.second);
      adj[id[e.first]].insert(id[e.second]);
    }

    inline bool
    contains_edge(const edge &e)
    {
      if (!contains_vertex(e.first) or !contains_vertex(e.second)) return false;
      auto id1 = id[e.first], id2 = id[e.second];
      return adj[id1].find(id2) != adj[id1].end();
    }
  };

}

