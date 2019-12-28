#pragma once

#include "base_graph.hpp"

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace simplex {

  // a directed, weighted graph
  template <typename Label_t=size_t, typename Weight_t = long long>
  struct AdjacencyList
    : public _GraphCRTP<AdjacencyList<Label_t, Weight_t>>
  {
    using _this_t = AdjacencyList<Label_t, Weight_t>;
    using _i_edge = typename GraphTraits<_this_t>::_i_edge;
    using _ie_hash = typename GraphTraits<_this_t>::_i_edge_hash;
    std::vector<std::unordered_set<_i_edge, _ie_hash>> adj;
    std::unordered_map<Label_t, size_t> id;
    std::unordered_map<size_t, Label_t> id_inv;
  public:
    using vertex = typename GraphTraits<_this_t>::vertex;
    using edge = typename GraphTraits<_this_t>::edge;
    using vertex_iter = typename GraphTraits<_this_t>::vertex_iter;
    using edge_iter = typename GraphTraits<_this_t>::edge_iter;
    using adjacent_iter = typename GraphTraits<_this_t>::adjacent_iter;
    using incident_iter = typename GraphTraits<_this_t>::incident_iter;

    AdjacencyList() {}

    void
    clear()
    {}

    iter_pair<vertex_iter>
    vertices() const
    { return iter_pair<vertex_iter>{vertex_iter{*this},
                                    vertex_iter{*this, true}}; }

    iter_pair<edge_iter>
    edges() const
    { return iter_pair<edge_iter>{edge_iter{*this},
                                  edge_iter{*this, true}}; }

    iter_pair<adjacent_iter>
    adjacent(const vertex &v) const
    { return iter_pair<adjacent_iter>{adjacent_iter{*this, v},
                                      adjacent_iter{*this, v, true}}; }

    iter_pair<incident_iter>
    incident(const vertex &v) const
    { return iter_pair<incident_iter>{incident_iter{*this, v},
                                      incident_iter{*this, v, true}}; }

    void
    add_vertex(const vertex &v)
    {
      if (id.find(v) != id.end()) return;
      auto k = id.size();
      id[v] = k;
      id_inv[k] = v;
      adj.emplace_back();
    }

    inline bool
    contains_vertex(const vertex &v) const
    { return id.find(v) != id.end(); }

    void
    add_edge(const edge &e)
    {
      if (contains_edge(e)) return;
      add_vertex(e.first);
      add_vertex(e.second);
      adj[id[e.first]].insert(_i_edge{id[e.second], e.weight});
    }

    inline bool
    contains_edge(const edge &e) const
    {
      if (!contains_vertex(e.first) or !contains_vertex(e.second)) 
        return false;
      auto id1 = id.at(e.first), id2 = id.at(e.second);
      return adj[id1].find(_i_edge{id2, e.weight}) != adj[id1].end();
    }

    size_t
    vertex_to_id(const vertex &v) const
    { return id.at(v); }

    const vertex&
    id_to_vertex(size_t id) const
    { return id_inv.at(id); }

    friend std::ostream&
    operator<<(std::ostream &o, const AdjacencyList<Label_t, Weight_t> &adj)
    {
      for (size_t i = 0; i < adj.adj.size(); ++i, o << "\n") {
        o << i << ": ";
        for (auto &x : adj.adj[i]) {
          o << x.to << ", ";
        }
      }
      return o;
    }
    
  private:
    friend class _GraphCRTP<AdjacencyList<Label_t, Weight_t>>;
    friend class GraphTraits<AdjacencyList<Label_t, Weight_t>>;
  };

  template <typename Label_t, typename Weight_t>
  struct GraphTraits<AdjacencyList<Label_t, Weight_t>> {
  private:
    using _that_t = AdjacencyList<Label_t, Weight_t>;
    struct _edge;
    struct _i_edge;
    class _vertex_iter;
    class _edge_iter;
    class _adjacent_iter;
    class _incident_iter;
  public:
    using vertex = Label_t;
    using edge = _edge;
    using vertex_iter = _vertex_iter;
    using edge_iter = _edge_iter;
    using adjacent_iter = _adjacent_iter;
    using incident_iter = _incident_iter;
  private:
    struct _edge {
      Label_t first, second;
      Weight_t weight;
    };

    struct _i_edge {
      size_t to;
      Weight_t weight;
      bool operator==(const _i_edge &other) const {
        return to == other.to and weight == other.weight;
      }
    };

    struct _i_edge_hash {
      auto operator()(const _i_edge &e) const {
        return std::hash<size_t>()(e.to + 280859 * e.weight);
      }
    };

    class _vertex_iter {
      typename std::unordered_map<Label_t, size_t>::const_iterator cur_it;
    public:
      // begin iter
      _vertex_iter(const _that_t &list): cur_it{list.id.begin()} {}

      // end iter
      _vertex_iter(const _that_t &list, bool): cur_it{list.id.end()} {}
      
      bool
      operator!=(const _vertex_iter &other) const
      { return cur_it != other.cur_it; }

      _vertex_iter&
      operator++()
      { ++cur_it; return *this; }

      const vertex&
      operator*() const
      { return cur_it->first; }
    };

    class _edge_iter {
      using adj_list_t = std::vector<std::unordered_set<_i_edge, _i_edge_hash>>;
      const _that_t *list;
      typename adj_list_t::const_iterator cur_u;
      typename std::unordered_set<_i_edge, _i_edge_hash>::const_iterator cur_v;
      
      bool
      _is_end() const
      {
        if (cur_u == list->adj.end()) return true;
        if (cur_u + 1 == list->adj.end() and cur_v == cur_u->end()) return true;
        return false;
      }
    public:
      explicit _edge_iter(const _that_t &list)
        : list{&list}, cur_u{list.adj.begin()}, cur_v{list.adj[0].begin()}
      { if (cur_v == cur_u->end()) this->operator++(); }

      explicit _edge_iter(const _that_t &list, bool)
        : list{&list}, cur_u{list.adj.end()}, cur_v{list.adj.back().begin()}
      {}

      _edge_iter(const _edge_iter &other)
        : list{other.list}, cur_u{other.cur_u}, cur_v{other.cur_v}
      {}

      _edge_iter&
      operator=(_edge_iter other) {
        using std::swap;
        swap(list, other.list);
        swap(cur_u, other.cur_u);
        swap(cur_v, other.cur_v);

        return *this;
      }

      bool
      operator!=(const _edge_iter &other) const
      {
        bool me = _is_end();
        bool you = other._is_end();
        if (me and you) return false;
        if ((me and !you) or (!me and you)) return true;
        if (list->adj != other.list->adj) return true;
        return cur_u != other.cur_u or cur_v != other.cur_v;
      }

      _edge_iter&
      operator++()
      {
        if (_is_end()) return *this;
        if (cur_v != cur_u->end() and ++cur_v != cur_u->end()) return *this;
        for (; cur_u != list->adj.end(); cur_v = (++cur_u)->begin()) {
          if (cur_v != cur_u->end()) return *this;
        }
        return *this;
      }

      _edge
      operator*()
      { return _edge{list->id_inv.at(cur_u - list->adj.begin()),
                     list->id_inv.at(cur_v->to),
                     cur_v->weight}; }
    };

    class _adjacent_iter {
      using iter_t = typename std::unordered_set<_i_edge, _i_edge_hash>::const_iterator;
      const _that_t *_list;
      iter_t it;
    public:
      _adjacent_iter(const _that_t &list, const vertex &v)
        : _list{&list}, it{list.adj[list.id.at(v)].begin()}
      {}

      _adjacent_iter(const _that_t &list, const vertex &v, bool)
        : _list{&list}, it{list.adj[list.id.at(v)].end()}
      {}
      
      bool
      operator!=(const _adjacent_iter &other) const
      { return it != other.it; }

      _adjacent_iter&
      operator++()
      { ++it; return *this; }

      vertex
      operator*()
      { return _list->id_inv.at(it->to); }
    };

    class _incident_iter {
      using iter_t = typename std::unordered_set<_i_edge, _i_edge_hash>::const_iterator;
      const _that_t *_list;
      iter_t it;
      vertex _from;
    public:
      _incident_iter(const _that_t &list, const vertex &v)
        : _list{&list}, it{list.adj[list.id.at(v)].begin()}, _from{v}
      {}

      _incident_iter(const _that_t &list, const vertex &v, bool)
        : _list{&list}, it{list.adj[list.id.at(v)].end()}, _from{v}
      {}
      
      bool
      operator!=(const _incident_iter &other) const
      { return it != other.it; }

      _incident_iter&
      operator++()
      { ++it; return *this; }

      edge
      operator*()
      { return edge{_from, _list->id_inv.at(it->to), it->weight}; }

    };

    friend class AdjacencyList<Label_t, Weight_t>;
  };
}

