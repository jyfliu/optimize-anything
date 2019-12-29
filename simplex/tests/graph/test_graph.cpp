#include <iostream>
#include "../../src/graph/graph.hpp"
#include "../../src/graph/graph_algos.hpp"

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
  ASSERT(graph.order() == 6);
  ASSERT(!graph.empty());
  ASSERT(graph.size() == 6);
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

  // TODO more silent asserts
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
auto test_bipartite_0() {
#define ASSERT(x) assert(x and "test_bipartite_0 ")
  using G = simplex::Undirected<simplex::AdjacencyList<std::string>>;
  G graph;
  graph.add_edge({"A", "a", 1});
  graph.add_edge({"B", "b", 2});
  graph.add_edge({"D", "c", 3});
  graph.add_edge({"C", "d", 4});
  graph.add_edge({"E", "a", 5});
  graph.add_edge({"A", "c", 6});
  graph.add_edge({"A", "e", 7});
  graph.add_edge({"B", "a", 8});
  graph.add_edge({"B", "b", 9});
  graph.add_edge({"C", "f", 10});
  graph.add_edge({"C", "a", 11});
  graph.add_edge({"D", "b", 12});
  graph.add_edge({"D", "c", 13});
  graph.add_edge({"F", "a", 14});
  graph.add_edge({"G", "a", 15});
  graph.add_edge({"H", "f", 16});
  graph.add_edge({"G", "g", 17});
  graph.add_edge({"H", "a", 18});
  graph.add_edge({"G", "b", 19});
  graph.add_edge({"C", "c", 20});
  graph.add_edge({"D", "a", 21});
  graph.add_edge({"D", "h", 22});
  simplex::Bipartition<G> bipartite(graph);
  ASSERT(simplex::num_connected_components(bipartite) == 1);
  ASSERT(bipartite.order() == 16);
  ASSERT(bipartite.size() == 22);
  ASSERT(bipartite.W().size() == 8);
  ASSERT(bipartite.J().size() == 8);
  bool b = 'A' <= (*bipartite.W().begin())[0] 
              and (*bipartite.W().begin())[0] <= 'Z';
  auto caps = b? bipartite.W() : bipartite.J();
  auto nots = b? bipartite.J() : bipartite.W();
  for (auto &&x : caps) ASSERT('A' <= x[0] and x[0] <= 'Z');
  for (auto &&y : bipartite.J()) ASSERT('a' <= y[0] and y[0] <= 'z');
#undef ASSERT
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
  test_bipartite_0();
  test_kruskals();

  std::cout << "All tests passed." << std::endl;
  return 0;
}

