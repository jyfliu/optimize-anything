#pragma once

#include "../Eigen/Dense"

#include <utility>
#include <vector>

namespace Simplex {
  class invalid_operation{};
  class size_mismatch{};
  /**
    * When result is
    *   solved => optimalSolution, optimalValue, certificate are defined
    *   unbounded => feasiblePoint, unboundedRay, certificate are defined
    *   infeasible => certificate is defined
    */
  class Result { 
    enum Status { solved, unbounded, infeasible };
    Status _status;
    Eigen::VectorXd _certificate;
    Eigen::VectorXd _vector;
    double _data;
    Result();
    Result(const Eigen::VectorXd &certificate);
    Result(const Eigen::VectorXd &certificate, const Eigen::VectorXd &vector);
    Result(const Eigen::VectorXd &certificate, const Eigen::VectorXd &vector,
           double data);
  public:
    static Result Solved(const Eigen::VectorXd &optimalSolution,
                         const Eigen::VectorXd &optimalCertificate,
                         double optimalValue);
    static Result Unbounded(const Eigen::VectorXd &initial,
                            const Eigen::VectorXd &ray);
    static Result Infeasible(const Eigen::VectorXd &infeasible);

    Status status() const noexcept;

    const Eigen::VectorXd &optimalSolution() const;
    double optimalValue() const;

    const Eigen::VectorXd &feasiblePoint() const;
    const Eigen::VectorXd &unboundedRay() const;

    const Eigen::VectorXd &certificate() const noexcept;
  };

  class LPProblem {
  public:
    const Eigen::MatrixXd A;
    const Eigen::VectorXd b;
    const Eigen::VectorXd c;
    LPProblem(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
              const Eigen::VectorXd &c);
  };

  Result solve(const LPProblem &problem);

  Result solveInteriorPoint(const LPProblem &problem);
  Result solveTwoPhaseSimplex(const LPProblem &problem);

}

