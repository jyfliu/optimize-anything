#include "eigen_utils.h"

#include "linear_programming_solver.h"


template <typename FieldType>
Simplex::Result<FieldType>::Result(
    const Simplex::Vector<FieldType> &certificate)
: _status{infeasible}, _certificate{certificate} {}

template <typename FieldType>
Simplex::Result<FieldType>::Result(
    const Simplex::Vector<FieldType> &certificate,
    const Simplex::Vector<FieldType> &vector)
: _status{unbounded}, _certificate{certificate}, _vector{vector} {}

template <typename FieldType>
Simplex::Result<FieldType>::Result(
    const Simplex::Vector<FieldType> &certificate,
    const Simplex::Vector<FieldType> &vector, FieldType data)
: _status{solved}, _certificate{certificate}, _vector{vector}, _data{data} {}

template <typename FieldType>
Simplex::Result<FieldType> Simplex::Result<FieldType>::Solved(
    const Simplex::Vector<FieldType> &optimalSolution,
    const Simplex::Vector<FieldType> &optimalCertificate,
    FieldType optimalValue)
{ return Simplex::Result<FieldType>(
      optimalCertificate, optimalSolution, optimalValue
  ); }

template <typename FieldType>
Simplex::Result<FieldType> Simplex::Result<FieldType>::Unbounded(
    const Simplex::Vector<FieldType> &initial,
    const Simplex::Vector<FieldType> &ray)
{ return Simplex::Result<FieldType>(ray, initial); }

template <typename FieldType>
Simplex::Result<FieldType> Simplex::Result<FieldType>::Infeasible(
    const Simplex::Vector<FieldType> &infeasible)
{ return Simplex::Result<FieldType>(infeasible); }

template <typename FieldType>
Simplex::Status Simplex::Result<FieldType>::status() const noexcept
{ return _status; }

template <typename FieldType>
const Simplex::Vector<FieldType>
&Simplex::Result<FieldType>::optimalSolution() const 
{ if (_status == solved) return _vector; else throw invalid_operation{}; }

template <typename FieldType>
FieldType Simplex::Result<FieldType>::optimalValue() const 
{ if (_status == solved) return _data; else throw invalid_operation{}; }

template <typename FieldType>
const Simplex::Vector<FieldType> &Simplex::Result<FieldType>::feasiblePoint() const 
{ if (_status == unbounded) return _vector; else throw invalid_operation{}; }

template <typename FieldType>
const Simplex::Vector<FieldType>
&Simplex::Result<FieldType>::unboundedRay() const
{ 
  if (_status == unbounded) return _certificate;
  else throw invalid_operation{};
}

template <typename FieldType>
const Simplex::Vector<FieldType>
&Simplex::Result<FieldType>::certificate() const noexcept
{ return _certificate; }

template <typename FieldType>
std::ostream &Simplex::operator<<(std::ostream &os,
                                  const Simplex::Result<FieldType> &r)
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



template <typename FieldType>
Simplex::LPProblem<FieldType>::LPProblem(const Matrix<FieldType> &A,
                                         const Vector<FieldType> &b,
                                         const Vector<FieldType> &c)
  : A{A}, b{b}, c{c}
{
  if (A.rows() != b.rows()) throw size_mismatch{};
  if (A.cols() != c.rows()) throw size_mismatch();
}

template <typename FieldType>
std::ostream &Simplex::operator<<(std::ostream &os,
                                  const Simplex::LPProblem<FieldType> &p)
{
  os << "LP: "; // TODO finish
  return os;
}


namespace {

  // simplex method with smallest subscripts rule / Bland's rule
  template <typename FieldType>
  Simplex::Result<FieldType> solveSimplexBland(
      const Simplex::LPProblem<FieldType> &problem,
      const Simplex::Vector<FieldType> &bfs,
      const Eigen::VectorXi &basis,
      int debugPrint = 0,
      FieldType epsilon = 1e-7)
  {
    // TODO add assertions to ensure input is clean
#define DPRINT(d) if (debugPrint >= d) std::cout
#define dout DPRINT(1)
    using namespace Eigen;
    using Simplex::Matrix;
    using Simplex::Vector;
    dout << "\nStarting new iteration" << std::endl;
    dout << "Basic feasible solution: x̅ = [" << bfs.transpose() << "]^T";
    dout << std::endl;
    dout << "Corresponding basis: B = {" << basis.transpose() << "}";
    dout << std::endl;
    Matrix<FieldType> A(problem.A);
    Vector<FieldType> b(problem.b);
    Vector<FieldType> c(problem.c);
    FieldType oValue = (c.transpose() * bfs)(0);
    dout << "Objective value: " << oValue << std::endl;
    const int n = A.cols();
    // N = [n] \ f
    VectorXi N(n - basis.rows());
    for (int i = 0, j = 0; i < n; ++i) {
      if (!Simplex::contains(basis, i)) N(j++) = i;
    }

    // definitions
    Matrix<FieldType> A_B = Simplex::subCols(A, basis);
    Matrix<FieldType> At_B = Simplex::subRows(A.transpose(), basis);
    Vector<FieldType> c_B = Simplex::subRows(c, basis);

    // solve A^T_B y=c_B
    dout << "Solving A^T_B y = c_B..." << std::endl;
    Vector<FieldType> y = At_B.colPivHouseholderQr().solve(c_B);
    dout << "Solution: y = [" << y.transpose() << "]^T" << std::endl;

    // compute c'_j = c_j - y^TA_j
    // let k be the first index st c'_k > 0
    // if no k exists then we have our optimal solution 
    // (every other direction is worse)
    int k = -1;
    for (int i = 0; i < N.rows(); ++i) {
      const int j = N(i);
      FieldType cp_j = c(j) - (y.transpose() * A.col(j))(0);
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
      return Simplex::Result<FieldType>::Solved(bfs, y, oValue);
    }
    dout << "Found entering index k = " << k << std::endl;
    dout << "Solving A_b d = A_k" << std::endl;
    VectorXi vec_k(1); vec_k << k;
    Vector<FieldType> A_k = Simplex::subCols(A, vec_k);
    Vector<FieldType> d = A_B.colPivHouseholderQr().solve(A_k);
    dout << "Solution: d = [" << d.transpose() << "]^T" << std::endl;

    // pick r with minimal x_r/d_r (=:xd)
    int r = -1;
    FieldType alpha = -1;
    for (int i = 0; i < d.rows(); ++i) {
      if (d(i) > epsilon) {
        FieldType xd = bfs(basis(i)) / d(i);
        if (r == -1 or xd < alpha - epsilon) {
          r = basis(i);
          alpha = xd;
        }
      }
    }
    Vector<FieldType> dp(n, 1);
    for (int j = 0, i = 0; j < n; ++j) {
      if (Simplex::contains(basis, j)) dp(j) = -d(i++);
      else if (j == k) dp(j) = 1;
      else dp(j) = 0;
    }
    if (r == -1) {
      dout << "No entry of d is positive, problem is unbounded." << std::endl;
      return Simplex::Result<FieldType>::Unbounded(bfs, dp);
    }
    dout << "The minimum ratio is α = " << alpha << std::endl;
    dout << "The leaving index is r = " << r << std::endl;
    Vector<FieldType> x = bfs + alpha * dp;
    VectorXi Bp(basis.rows());
    int i = 0;
    for (int j = 0; j < n; ++j) {
      if (j != r and (Simplex::contains(basis, j) or j == k)) {
        Bp(i++) = j;
      }
    }
    assert(i == basis.rows() and "something went wrong, basis creation ");
    return solveSimplexBland(problem, x, Bp, debugPrint, epsilon);
#undef DPRINT
#undef dout
  }
}

template <typename FieldType>
Simplex::Result<FieldType> Simplex::solve(const Simplex::LPProblem<FieldType> &problem)
{ // TODO: use different solvers
  return Simplex::solveTwoPhaseSimplex(problem);
}

template <typename FieldType>
Simplex::Result<FieldType> Simplex::solveTwoPhaseSimplex(
    const Simplex::LPProblem<FieldType> &problem,
    int debugPrint, FieldType epsilon)
{ // TODO write a version which works for rationals. Probably can not use qr_decomp
#define DPRINT(d) if (debugPrint >= d) std::cout
#define dout DPRINT(1)
  using namespace Eigen;
  // Ax = b, x >= 0
  // A = QRP^T
  // QRP^Tx = b, x >= 0
  // PHASE 1 problem: max [0 -1]^Ty st [R 1]y = Q^-1b, y >= 0
  Matrix<FieldType> A = problem.A;
  Vector<FieldType> b = problem.b;
  Vector<FieldType> c = problem.c;

  auto qr_decomp = A.colPivHouseholderQr();
  int rank = qr_decomp.rank();
  Matrix<FieldType> Q(qr_decomp.matrixQ());
  // ensure that Ap has full row rank
  Matrix<FieldType> R(qr_decomp.matrixQR()
      .topRows(rank)
      .template triangularView<Upper>());
  Matrix<FieldType> P(qr_decomp.colsPermutation());
  dout << Q << std::endl << R << std::endl << P << std::endl;

  Matrix<FieldType> Ap(R.rows(), R.cols() + rank);
  for (int i = 0; i < R.cols(); ++i) { Ap.col(i) = R.col(i); }
  Vector<FieldType> bp = (Q.inverse() * b).head(rank);
  dout << "b\n" << bp << std::endl;
  // ensure b has only non-negative entries
  for (int i = 0; i < rank; ++i) {
    if (bp(i) < 0) {
      bp(i) *= -1;
      Ap.row(i) *= -1;
    }
  }
  for (int i = 0; i < rank; ++i) {
    Ap.col(i + R.cols()) = Vector<FieldType>::Unit(rank, i);
  }
  Vector<FieldType> c_aux(R.cols() + rank);
  for (int i = 0; i < R.cols(); ++i) { c_aux(i) = 0; }
  for (int i = R.cols(); i < R.cols() + rank; ++i) { c_aux(i) = -1; }
  VectorXi basis_aux(rank);
  for (int i = 0; i < rank; ++i) { basis_aux(i) = i + R.cols(); }
  Vector<FieldType> bfs_aux(R.cols() + rank);
  for (int i = 0; i < R.cols(); ++i) { bfs_aux(i) = 0; }
  for (int i = 0; i < rank; ++i) { bfs_aux(i + R.cols()) = bp(i); }

  dout << "A:\n" << Ap << "\nb\n" << bp << "\nc\n" << c_aux << std::endl;
  Simplex::LPProblem<FieldType> problem_aux(Ap, bp, c_aux);
  Simplex::Result<FieldType> phase1Result = solveSimplexBland(
    problem_aux, bfs_aux, basis_aux, debugPrint, epsilon
  );
  dout << phase1Result << std::endl;
  if (phase1Result.optimalValue() < - epsilon) {
    //TODO transform the certificate to obtain a certificate for the original problem
    return Simplex::Result<FieldType>::Infeasible(phase1Result.certificate());
  }
  // PHASE 2 problem: max c^TPy st Ry = Q^-1b, y >= 0
  Vector<FieldType> bfs_main = phase1Result.optimalSolution().head(R.cols());
  VectorXi basis_main(rank);
  for (int i = 0, j = 0; i < R.cols(); ++i) {
    if (bfs_main(i) > epsilon) basis_main(j++) = i;
  }
  dout << bfs_main << std::endl << basis_main << std::endl;
  Vector<FieldType> c_main = (c.transpose() * P).transpose();
  dout << R << "\n" << bp << std::endl;
  Simplex::LPProblem<FieldType> problem_main(R, bp, c_main);
  return solveSimplexBland(
    problem_main, bfs_main, basis_main, debugPrint, epsilon
  );
#undef DPRINT
#undef dout
}

