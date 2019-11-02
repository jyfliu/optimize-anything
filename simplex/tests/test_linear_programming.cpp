#include <iostream>
#include <cmath>

#include "../src/utils/rational.h"
#include "../src/linear_programming.cpp"


constexpr double epsilon = 1e-7;

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
  Simplex::LPProblem<double> p(A, b, c);
  Simplex::Result<double> res = solveSimplexBland(p, bfs, basis, 0);
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

void test_solveSimplexBlandRational() {
#define ASSERT(x) assert(x and "test_solveSimplexBlandRational ")
  typedef Simplex::rational<int> Q;
  Eigen::Matrix<Q, -1, -1> A(4, 7);
  A << 1, 2, 2, 2, 2, 2, 2,
       1, 1, 1, 2, 1, 1, 1,
       1, 1, 2, 1, 4, -2, 0,
       1, 2, 1, 1, -1, 4, 0;
  Eigen::Matrix<Q, -1, 1> b(4);
  b << 7, 5, 5, 5;
  Eigen::Matrix<Q, -1, 1> c(7);
  c << 1, 2, 2, 2, 6, 4, -4;
  Eigen::VectorXi basis(4);
  basis << 0, 1, 2, 3;
  Eigen::Matrix<Q, -1, 1> bfs(7);
  bfs << 1, 1, 1, 1, 0, 0, 0;
  Simplex::LPProblem<Q> p(A, b, c);
  Simplex::Result<Q> res = solveSimplexBland<Q>(p, bfs, basis, 0, Q(0),
      Simplex::FullPivLUSolver<Q>);
  //std::cout << res << std::endl;
  ASSERT(res.status() == Simplex::Status::solved);
  Q opt = Q(145, 11);
  //Eigen::VectorXd sol(7);
  //sol << 1.9090909090, 0, 0, 0.5454545454, 1.0909090909, 0.9090909090, 0;
  //Eigen::VectorXd cer(4);
  //cer << 4.0909090909, -3.0909090909, 0.1818181818, -0.1818181818;
  ASSERT(res.optimalValue() == opt);
  //for (int i = 0; i < 7; ++i)
    //ASSERT(std::abs(res.optimalSolution()(i) - sol(i)) < epsilon);
  //for (int i = 0; i < 4; ++i)
    //ASSERT(std::abs(res.certificate()(i) - cer(i)) < epsilon);
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
  Simplex::LPProblem<double> p(A, b, c);
  Simplex::Result<double> res = solveTwoPhaseSimplexQR(p, 0);
  //std::cout << res << std::endl;
  ASSERT(res.status() == Simplex::Status::solved);
  double opt = 13.1818181818;
  ASSERT(std::abs(res.optimalValue() - opt) < epsilon);
#undef ASSERT
}

void test_solveTwoPhaseSimplexRational() {
#define ASSERT(x) assert(x and "test_solveTwoPhaseSimplexRational ")
  typedef Simplex::rational<int> Q;

  Eigen::Matrix<Q, -1,-1> A(7, 7);
  A << 1, 2, 2, 2, 2, 2, 2,
       1, 1, 1, 2, 1, 1, 1,
       1, 1, 2, 1, 4, -2, 0,
       1, 2, 1, 1, -1, 4, 0,
       2, 2, 3, 3, 5, -1, 1,
       2, 3, 3, 4, 3, 3, 3,
       0, 0, 1, -1, 3, -3, -1;

  Eigen::Matrix<Q, -1,1> b(7);
  b << 7, 5, 5, 5, 10, 12, 0;
  Eigen::Matrix<Q, -1, 1> c(7);
  c << 1, 2, 2, 2, 6, 4, -4;
  Simplex::LPProblem<Q> p(A, b, c);
  Simplex::Result<Q> res = solveTwoPhaseSimplexLU(p, 0, Q(0));
  //std::cout << res << std::endl;
  ASSERT(res.status() == Simplex::Status::solved);
  Q opt = Q(145, 11);
  ASSERT(res.optimalValue() == opt);
#undef ASSERT
}


int main() {
  test_solveSimplexBland();
  test_solveSimplexBlandRational();
  test_solveTwoPhaseSimplex();
  test_solveTwoPhaseSimplexRational();
  std::cout << "All tests passed." << std::endl;
  return 0;
}

