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

Simplex::Status Simplex::Result::status() const noexcept
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

std::ostream &operator<<(std::ostream &os, const Simplex::Result &r)
{
  os << "Result: ";
  switch (r.status()) {
    case Simplex::Status::solved:
      os << "solved" << std::endl;
      os << "optimal solution:\n" << r.optimalSolution() << std::endl;
      os << "optimal value: " << r.optimalValue() << std::endl;
      os << "optimality certificate:\n" << r.certificate() << std::endl;
      break;
    case Simplex::Status::unbounded:
      os << "unbounded" << std::endl;
      os << "feasible solution:\n" << r.optimalSolution() << std::endl;
      os << "ray certificate :\n" << r.certificate() << std::endl;
      break;
    case Simplex::Status::infeasible:
      os << "infeasible\n" << std::endl;
      os << "infeasibility certificate :\n" << r.certificate() << std::endl;
      break;
  }
  return os;
}



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
    for (int i = 0; i < basis.rows(); ++i) {
      toBeReturned.row(i) = mat.row(basis(i));
    }
    return toBeReturned;
  }

  Eigen::MatrixXd subCols(const Eigen::MatrixXd &mat,
                          const Eigen::VectorXi &basis)
  {
    Eigen::MatrixXd toBeReturned(mat.rows(), basis.rows());
    for (int i = 0; i < basis.rows(); ++i) {
      toBeReturned.col(i) = mat.col(basis(i));
    }
    return toBeReturned;
  }

  bool contains(const Eigen::MatrixXi& Z, int d)
  {
    for (int i = 0; i < Z.rows(); ++i)
      for (int j = 0; j < Z.cols(); ++j)
        if (Z(i, j) == d) return true;
    return false;
  }

  // simplex method with smallest subscripts rule / Bland's rule
  Simplex::Result solveSimplexBland(const Simplex::LPProblem &problem,
                                    const Eigen::VectorXd &bfs,
                                    const Eigen::VectorXi &basis,
                                    int debugPrint = 0,
                                    int epsilon = 1e-7)
  {
#define DPRINT(d) if (debugPrint >= d) std::cout
#define dout DPRINT(1)
    dout << "Basic feasible solution: x̅ = [" << bfs.transpose() << "]^T";
    dout << std::endl;
    dout << "Corresponding basis: B = {" << basis.transpose() << "}";
    dout << std::endl;
    using namespace Eigen;
    MatrixXd A(problem.A);
    VectorXd b(problem.b);
    VectorXd c(problem.c);
    const int n = A.cols();
    // N = [n] \ f
    VectorXi N(n - basis.rows());
    for (int i = 0, j = 0; i < n; ++i) {
      if (!contains(basis, i)) N(j++) = i;
    }

    // definitions
    MatrixXd A_B = subCols(A, basis);
    MatrixXd At_B = subRows(A.transpose(), basis);
    VectorXd c_B = subRows(c, basis);

    // solve A^T_B y=c_B
    dout << "Solving A^T_B y = c_B..." << std::endl;
    VectorXd y = At_B.colPivHouseholderQr().solve(c_B);
    dout << "Solution: y = [" << y.transpose() << "]^T\n" << std::endl;

    // compute c'_j = c_j - y^TA_j
    // let k be the first index st c'_k > 0
    // if no k exists then we have our optimal solution 
    // (every other direction is worse)
    int k = -1;
    for (int i = 0; i < N.rows(); ++i) {
      int j = N(i);
      int cp_j = c(j) - (y.transpose() * A.col(j))(0);
      dout << "Computing c̅_" << j << " = c_" << j << " - y^T A_" << j;
      dout << std::endl;
      dout << "Solution: c̅_" << j << " = " <<  cp_j << std::endl;
      if (cp_j > epsilon) {
        k = j;
        break;
      }
    }
    // we found optimal solution
    if (k == -1) {
      dout << "No entering value exists, x̅ is optimal." << std::endl;
      return Simplex::Result::Solved(bfs, y, (c.transpose() * bfs)(0));
    }
    dout << "Found entering indice k = " << k << std::endl;
    dout << "Solving A_b d = A_k" << std::endl;
    VectorXi vec_k(1); vec_k << k;
    VectorXd A_k = subCols(A, vec_k);
    VectorXd d = A_B.colPivHouseholderQr().solve(A_k);
    dout << "Solution: d = [" << d.transpose() << "]^T" << std::endl;

    // pick r with minimal x_r/d_r (=:xd)
    int r = -1;
    double alpha = -1;
    for (int i = 0; i < d.rows(); ++i) {
      if (d(i) > epsilon) {
        double xd = bfs(basis(i)) / d(i);
        if (r == -1 or xd < alpha - epsilon) {
          r = basis(i);
          alpha = xd;
        }
      }
    }
    VectorXd dp(n, 1);
    for (int j = 0, i = 0; j < n; ++j) {
      if (contains(basis, j)) dp(j) = -d(i++);
      else if (j == k) dp(j) = 1;
      else dp(j) = 0;
    }
    if (r == -1) {
      dout << "No entry of d is positive, problem is unbounded." << std::endl;
      return Simplex::Result::Unbounded(bfs, dp);
    }
    dout << "The minimum ratio is α = " << alpha << std::endl;
    dout << "The leaving index is r = " << r << std::endl;
    VectorXd x = bfs + alpha * dp;
    VectorXi Bp(basis.rows());
    int i = 0;
    for (int j = 0; j < n; ++j) {
      if (j != r and (contains(basis, j) or j == k)) {
        Bp(i++) = j;
      }
    }
    assert(i == basis.rows() and "something went wrong, basis creation ");
    return solveSimplexBland(problem, x, Bp, debugPrint, epsilon);
#undef DPRINT
#undef dout
  }
}

Simplex::Result Simplex::solve(const Simplex::LPProblem &problem)
{ // TODO: use different solvers
  return Simplex::solveTwoPhaseSimplex(problem);
}

Simplex::Result Simplex::solveTwoPhaseSimplex(const Simplex::LPProblem &problem)
{ throw 0; }

