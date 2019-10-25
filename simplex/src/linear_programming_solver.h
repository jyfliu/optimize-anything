#pragma once

#include "../Eigen/Dense"

#include <utility>

namespace Simplex {
  class Result { 
    enum Status { solved, unbounded, infeasible };
    Status status;
    std::pair<Eigen::VectorXd, double> solution;
    std::pair<Eigen::VectorXd, Eigen::VectorXd> unboundedCert;
    Eigen::VectorXd infeasibleCert;
    Result();
  public:
    static Result Solved(const Eigen::VectorXd &optimalSolution, double optimalValue);
    static Result Unbounded(const Eigen::VectorXd &initial, const Eigen::VectorXd &ray);
    static Result Infeasible(const Eigen::VectorXd &infeasible);
  };

  class LPProblem {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd c;

  public:
    LPProblem(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
              const Eigen::VectorXd &c);
  };

  Result solve(const LPProblem &problem);

  Result solveInteriorPoint(const LPProblem &problem);
  Result solveSimplex(const LPProblem &problem);

}

