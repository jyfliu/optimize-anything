#pragma once

#include <unordered_map>
#include <vector>

namespace simplex {

  template <typename T>
  class disjoint_set {
    std::unordered_map<T, size_t> id;
    //std::unordered_map<size_t, T> id_inv;
    std::vector<size_t> p;
    std::vector<size_t> r;
  public:
    disjoint_set() {}

    disjoint_set(size_t size)
    { p.reserve(size); r.reserve(size); }

    void
    make_set(const T &t)
    {
      if (id.find(t) != id.end()) return;
      auto k = id.size();
      id[t] = k;
      //id_inv[k] = t;
      p.push_back(k);
      r.push_back(0);
    }

    size_t
    find(const T &t)
    {
      auto &it = id.find(t);
      if (it == id.end()) return -1;
      return _find(*it);
    }

    inline size_t
    end()
    { return -1; }

    void
    union_set(const T &u, const T &v)
    {
      auto &it1 = id.find(u), &it2 = id.find(v);
      if (it1 == id.end() or it2 == id.end()) return;
      auto x = _find(*it1), y = _find(*it2);
      if (r[x] > r[y]) p[y] = x;
      else if (r[y] > r[x]) p[x] = y;
      else {
        p[x] = y;
        ++r[y];
      }
    }
  private:
    inline size_t
    _find(size_t x)
    {
      if (p[x] != x) p[x] = _find(p[x]);
      return p[x];
    }
  };
}

