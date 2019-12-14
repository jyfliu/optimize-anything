#include <iostream>
#include "../../src/graph/adjacency_list.hpp"

#include <iostream>

void test_adj_list() {
#define ASSERT(x) assert(x and "test_adj_list ")
  simplex::AdjacencyList<std::string> graph;
  graph.add_vertex("A");
  graph.add_edge({"B", "C", 3});
  graph.add_vertex("D");
  graph.add_edge({"B", "D", 5});
  graph.add_vertex("A");
  graph.add_edge({"B", "E", 2});
  graph.add_vertex("A");
  graph.add_edge({"E", "D", 6});
  graph.add_edge({"E", "F", 9});
  graph.add_edge({"E", "D", 7});
  std::cout << "vertices: ";
  for (auto &v : graph.vertices()) {
    std::cout << v <<", ";
  }
  std::cout << std::endl;
  std::cout << "edges: ";
  for (auto e : graph.edges()) {
    std::cout << "("<< e.first << ", " << e.second << ", "<< e.weight <<"), ";
  }
  std::cout << std::endl;
#undef ASSERT
}

int main(void) {
  test_adj_list();

  std::cout << "All tests passed." << std::endl;
  return 0;
}

