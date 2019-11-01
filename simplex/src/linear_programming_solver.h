#pragma once

#include "eigen_utils.h"
#include "../Eigen/Dense"

#include <utility>
#include <vector>
#include <iostream>
#include <functional>

namespace Simplex {
  class invalid_operation{};
  class size_mismatch{};

  // typedef a vector of field members
  template <typename FieldType>
  using Vector = Eigen::Matrix<FieldType, Eigen::Dynamic, 1>;

  template <typename FieldType>
  using Matrix = Eigen::Matrix<FieldType, Eigen::Dynamic, Eigen::Dynamic>;

  enum Status { solved, unbounded, infeasible };
  /**
    * When result is
    *   solved => optimalSolution, optimalValue, certificate are defined
    *   unbounded => feasiblePoint, unboundedRay, certificate are defined
    *   infeasible => certificate is defined
    */
  template <typename FieldType>
  class Result {
    Status _status;
    Vector<FieldType> _certificate;
    Vector<FieldType> _vector;
    FieldType _data;
    Result(const Vector<FieldType> &certificate);
    Result(const Vector<FieldType> &certificate,
           const Vector<FieldType> &vector);
    Result(const Vector<FieldType> &certificate,
           const Vector<FieldType> &vector,
           FieldType data);
  public:
    static Result Solved(const Vector<FieldType> &optimalSolution,
                         const Vector<FieldType> &optimalCertificate,
                         FieldType optimalValue);
    static Result Unbounded(const Vector<FieldType> &initial,
                            const Vector<FieldType> &ray);
    static Result Infeasible(const Vector<FieldType> &infeasible);

    Status status() const noexcept;

    const Vector<FieldType> &optimalSolution() const;
    FieldType optimalValue() const;

    const Vector<FieldType> &feasiblePoint() const;
    const Vector<FieldType> &unboundedRay() const;

    const Vector<FieldType> &certificate() const noexcept;
  };

  // an LP in SEF. TODO refactor and support LPs in other forms
  template <typename FieldType>
  class LPProblem {
  public:
    const Matrix<FieldType> A;
    const Vector<FieldType> b;
    const Vector<FieldType> c;
    LPProblem(const Matrix<FieldType> &A, const Vector<FieldType> &b,
              const Vector<FieldType> &c);
  };

  template <typename FieldType>
  Result<FieldType> solve(const LPProblem<FieldType> &problem);

  template <typename FieldType>
  Result<FieldType> solveInteriorPoint(const LPProblem<FieldType> &problem);

  template <typename FieldType>
  Result<FieldType> solveTwoPhaseSimplexQR(const LPProblem<FieldType> &problem, 
                                           int debugPrint=0,
                                           FieldType epsilon=1e-7);

  template <typename FieldType>
  Result<FieldType> solveTwoPhaseSimplexLU(const LPProblem<FieldType> &problem, 
                                           int debugPrint=0,
                                           FieldType epsilon=1e-7);

  template <typename FieldType>
  std::ostream &operator<<(std::ostream &os,
                           const Simplex::Result<FieldType> &r);

  template <typename FieldType>
  std::ostream &operator<<(std::ostream &os,
                           const Simplex::LPProblem<FieldType> &p);

}

