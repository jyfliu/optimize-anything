#include "../src/linear_programming_solver.cpp"

#include <iostream>
#include <cmath>

constexpr double epsilon = 1e-7;

struct test_failed{ std::string msg; };

void test_subRows() {
#define ASSERT(x) assert(x and "test_subRows ")
  Eigen::MatrixXd mat(5, 5);
  mat << 1, 2, 3, 4, 5,
         6, 7, 8, 9, 10,
         11, 12, 13, 14, 15,
         16, 17, 18, 19, 20,
         21, 22, 23, 24, 25;
  Eigen::VectorXi vec(3);
  vec << 0, 2, 3;
  Eigen::MatrixXd res(3, 5);
  res << 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
  ASSERT(subRows(mat, vec) == res);
#undef ASSERT
}

void test_solveSimplexBland() {
#define ASSERT(x) assert(x and "test_solveSimplexBland ")
  Eigen::MatrixXd A(4, 7);
  A << 1, 2, 2, 2, 2, 2, 2,
       1, 1, 1, 2, 1, 1, 1,
       1, 1, 2, 1, 4, -2, 0,
       1, 2, 1, 1, -1, 4, 0;
  Eigen::VectorXd b(4);
  b << 7, 5, 5, 5;
  Eigen::VectorXd c(7);
  c << 1, 2, 2, 2, 6, 4, -4;
  Eigen::VectorXi basis(4);
  basis << 0, 1, 2, 3;
  Eigen::VectorXd bfs(7);
  bfs << 1, 1, 1, 1, 0, 0, 0;
  Simplex::LPProblem p(A, b, c);
  Simplex::Result res = solveSimplexBland(p, bfs, basis, 0);
  //std::cout << res << std::endl;
  ASSERT(res.status() == Simplex::Status::solved);
  double opt = 13.1818181818;
  Eigen::VectorXd sol(7);
  sol << 1.9090909090, 0, 0, 0.5454545454, 1.0909090909, 0.9090909090, 0;
  Eigen::VectorXd cer(4);
  cer << 4.0909090909, -3.0909090909, 0.1818181818, -0.1818181818;
  ASSERT(std::abs(res.optimalValue() - opt));
  for (int i = 0; i < 7; ++i)
    ASSERT(std::abs(res.optimalSolution()(i) - sol(i)) < epsilon);
  for (int i = 0; i < 4; ++i)
    ASSERT(std::abs(res.certificate()(i) - cer(i)) < epsilon);
#undef ASSERT
}

void test_solveTwoPhaseSimplex() {
#define ASSERT(x) assert(x and "test_solveTwoPhaseSimplex ")
  Eigen::MatrixXd A(7, 7);
  A << 1, 2, 2, 2, 2, 2, 2,
       1, 1, 1, 2, 1, 1, 1,
       1, 1, 2, 1, 4, -2, 0,
       1, 2, 1, 1, -1, 4, 0,
       2, 2, 3, 3, 5, -1, 1,
       2, 3, 3, 4, 3, 3, 3,
       0, 0, 1, -1, 3, -3, -1;

  Eigen::VectorXd b(7);
  b << 7, 5, 5, 5, 10, 12, 0;
  Eigen::VectorXd c(7);
  c << 1, 2, 2, 2, 6, 4, -4;
  Simplex::LPProblem p(A, b, c);
  Simplex::Result res = solveTwoPhaseSimplex(p, 0);
  //std::cout << res << std::endl;
  ASSERT(res.status() == Simplex::Status::solved);
  double opt = 13.1818181818;
  ASSERT(std::abs(res.optimalValue() - opt));
#undef ASSERT
}


int main() {
  test_subRows();
  test_solveSimplexBland();
  test_solveTwoPhaseSimplex();
  std::cout << "All tests passed." << std::endl;
  return 0;
}

