#include <iostream>
#include <cassert>

#include "../../src/graph/depth_first_search.hpp"
#include "../../src/graph/graph.hpp"

void test_basic_0() {
#define ASSERT(x) assert(x and "test_basic_0 ")
  simplex::Graph<std::string> g;
  g.add_edge({"A", "B", 2});
  g.add_edge({"A", "D", 2});
  g.add_edge({"A", "E", 2});
  g.add_edge({"B", "C", 2});
  g.add_edge({"F", "C", 2});
  g.add_edge({"H", "D", 2});
  g.add_edge({"I", "F", 2});
  g.add_edge({"G", "A", 2});
  
  struct Visitor {
    size_t cnt = 0;
    void discover_vertex(const std::string&, const simplex::Graph<std::string>&)
    { ++cnt; }
  } visitor;

  auto res = simplex::depth_first_search(g, "A", visitor);
  ASSERT(res.cnt == g.order());
#undef ASSERT
}

int main() {
  test_basic_0();
  std::cout << "All tests passed." << std::endl;
  return 0;
}

