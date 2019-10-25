#include "linear_programming_solver.h"


Simplex::Result::Result(const Eigen::VectorXd &certificate):
  _status{infeasible}, _certificate{certificate} {}

Simplex::Result::Result(const Eigen::VectorXd &certificate,
                        const Eigen::VectorXd &vector):
  _status{unbounded}, _certificate{certificate}, _vector{vector} {}

Simplex::Result::Result(const Eigen::VectorXd &certificate,
                        const Eigen::VectorXd &vector, double data):
  _status{solved}, _certificate{certificate}, _vector{vector}, _data{data} {}

Simplex::Result Simplex::Result::Solved(
    const Eigen::VectorXd &optimalSolution,
    const Eigen::VectorXd &optimalCertificate,
    double optimalValue)
{ return Simplex::Result(optimalCertificate, optimalSolution, optimalValue); }

Simplex::Result Simplex::Result::Unbounded(
    const Eigen::VectorXd &initial, const Eigen::VectorXd &ray)
{ return Simplex::Result(ray, initial); }

Simplex::Result Simplex::Result::Infeasible(const Eigen::VectorXd &infeasible)
{ return Simplex::Result(infeasible); }

Simplex::Result::Status Simplex::Result::status() const noexcept
{ return _status; }

const Eigen::VectorXd &Simplex::Result::optimalSolution() const 
{ if (_status == solved) return _vector; else throw invalid_operation{}; }

double Simplex::Result::optimalValue() const 
{ if (_status == solved) return _data; else throw invalid_operation{}; }

const Eigen::VectorXd &Simplex::Result::feasiblePoint() const 
{ if (_status == unbounded) return _vector; else throw invalid_operation{}; }

const Eigen::VectorXd &Simplex::Result::unboundedRay() const
{ 
  if (_status == unbounded) return _certificate;
  else throw invalid_operation();
}

const Eigen::VectorXd &Simplex::Result::certificate() const noexcept
{ return _certificate; }



Simplex::LPProblem::LPProblem(const Eigen::MatrixXd &A,
                              const Eigen::VectorXd &b,
                              const Eigen::VectorXd &c): A{A}, b{b}, c{c}
{
  if (A.rows() != b.rows()) throw size_mismatch{};
  if (A.cols() != c.rows()) throw size_mismatch();
}



namespace {
  
  Eigen::MatrixXd subRows(const Eigen::MatrixXd &mat,
                          const Eigen::VectorXi &basis)
  {
    Eigen::MatrixXd toBeReturned(basis.rows(), mat.cols());
    for (int i=0; i<basis.rows(); ++i) {
      toBeReturned.row(i) = mat.row(basis(i));
    }
    return toBeReturned;
  }

  // simplex method with smallest subscripts rule / Bland's rule
  Simplex::Result solveSimplexBland(const Simplex::LPProblem &problem,
                                    const Eigen::VectorXd &bfs,
                                    const Eigen::VectorXd &basis)
  {
    using namespace Eigen;
    MatrixXd A(problem.A);
    VectorXd b(problem.b);
    VectorXd c(problem.c);
    const int n = A.rows();
    const int m = A.cols();
    // N = [n] \ f
    VectorXd N = MatrixXd::Zero(n - basis.rows(), 1);
    for (int i=0, j=0; i<n; i++) {
      if ((basis.array() != i).any()) N(j++) = i;
    }

    // definitions
    MatrixXd A_B = subRows(A, basis);
    MatrixXd At_B = subRows(A.transpose(), basis);
    VectorXd c_B = subRows(c, basis);

    // solve A^T_B y=c_B
    VectorXd y = At_B.colPivHouseholderQr().solve(c);

    // TODO
  }
}

Simplex::Result Simplex::solve(const Simplex::LPProblem &problem)
{ // TODO: use different solvers
  return Simplex::solveTwoPhaseSimplex(problem);
}

Simplex::Result Simplex::solveTwoPhaseSimplex(const Simplex::LPProblem &problem)
{ throw 0; }

