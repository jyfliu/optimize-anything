#include <iostream>
#include "../../src/graph/adjacency_list.hpp"
#include "../../src/graph/undirected_graph.hpp"
#include "../../src/graph/minimum_spanning_tree.hpp"

#include <iostream>
#include <cassert>

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
  // TODO have silent tests
  //std::cout << "vertices: ";
  //for (auto &v : graph.vertices()) {
    //std::cout << v <<", ";
  //}
  //std::cout << std::endl;
  //std::cout << "edges: ";
  //for (auto e : graph.edges()) {
    //std::cout << "("<< e.first << ", " << e.second << ", "<< e.weight <<"), ";
  //}
  //std::cout << std::endl;
  //std::cout << "adj: ";
  //for (auto v : graph.vertices()) {
    //std::cout << v << " ";
    //for (auto v2 : graph.adjacent(v)) {
      //std::cout << v2 << " ";
    //}
    //std::cout << std::endl;
  //}
  //std::cout << std::endl;
  //std::cout << "inc: ";
  //for (auto v : graph.vertices()) {
    //std::cout << v << " ";
    //for (auto e : graph.incident(v)) {
      //std::cout << e.second << " ";
    //}
    //std::cout << std::endl;
  //}
#undef ASSERT
}

auto test_undirected() {
  // also test kruskals
#define ASSERT(x) assert(x and "test_undirected ")
  simplex::Undirected<simplex::AdjacencyList<std::string>> g;
  g.add_edge({"A", "B", 2});
  g.add_edge({"B", "C", 3});
  g.add_edge({"C", "D", 1});
  g.add_edge({"D", "E", 15});
  g.add_edge({"E", "A", 17});
  g.add_edge({"A", "a", 8});
  g.add_edge({"B", "b", 11});
  g.add_edge({"C", "c", 9});
  g.add_edge({"D", "d", 16});
  g.add_edge({"E", "e", 7});
  g.add_edge({"a", "c", 3});
  g.add_edge({"b", "d", 19});
  g.add_edge({"c", "e", 2});
  g.add_edge({"d", "a", 6});
  g.add_edge({"e", "b", 4});

  // TODO silent asserts
  //for (auto v : g.vertices()) {
    //std::cout << v << ": ";
    //for (auto v2 : g.adjacent(v)) {
      //std::cout << v2 << ", ";
    //}
    //std::cout << std::endl;
  //}
  //for (auto e : g.edges()) {
    //std::cout << e.first << ", "<< e.second << ", " << e.weight << std::endl;
  //}
#undef ASSERT
  return g;
}
void test_kruskals() {
  // also test kruskals
#define ASSERT(x) assert(x and "test_kruskals ")
  auto g = test_undirected();
  auto h = kruskals(g);
  long long s = 0;
  long long c = 0;
  for (auto &&e : h.edges()) {
    s += e.weight;
    ++c;
  }
  ASSERT(s == 36);
  ASSERT(c == 9);
#undef ASSERT
}
int main(void) {
  test_adj_list();
  test_undirected();
  test_kruskals();

  std::cout << "All tests passed." << std::endl;
  return 0;
}

