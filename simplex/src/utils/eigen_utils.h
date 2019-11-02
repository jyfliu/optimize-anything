#include "../../Eigen/Dense"

namespace Simplex {

  // typedef a vector of field members
  template <typename FieldType>
  using Vector = Eigen::Matrix<FieldType, Eigen::Dynamic, 1>;

  template <typename FieldType>
  using Matrix = Eigen::Matrix<FieldType, Eigen::Dynamic, Eigen::Dynamic>;

  template <typename FieldType>
  using Solver = std::function<Simplex::Vector<FieldType>
    (Simplex::Matrix<FieldType>, Simplex::Vector<FieldType>)>;

  // what the f*ck eigen??
  template <typename Derived>
  Eigen::Matrix<
    typename Eigen::internal::traits<Derived>::Scalar,
    Eigen::internal::traits<Derived>::RowsAtCompileTime, 
    Eigen::internal::traits<Derived>::ColsAtCompileTime>
  subRows(const Derived &mat,
          const Eigen::VectorXi &basis)
  {
    Eigen::Matrix<typename Eigen::internal::traits<Derived>::Scalar,
     Eigen::internal::traits<Derived>::RowsAtCompileTime, 
     Eigen::internal::traits<Derived>::ColsAtCompileTime>
    toBeReturned(basis.rows(), mat.cols());

    for (int i = 0; i < basis.rows(); ++i) {
      toBeReturned.row(i) = mat.row(basis(i));
    }
    return toBeReturned;
  }

  template <typename Derived>
  Eigen::Matrix<
    typename Eigen::internal::traits<Derived>::Scalar,
    Eigen::internal::traits<Derived>::RowsAtCompileTime, 
    Eigen::internal::traits<Derived>::ColsAtCompileTime>
  subCols(const Derived &mat,
          const Eigen::VectorXi &basis)
  {
    Eigen::Matrix<typename Eigen::internal::traits<Derived>::Scalar,
     Eigen::internal::traits<Derived>::RowsAtCompileTime, 
     Eigen::internal::traits<Derived>::ColsAtCompileTime>
    toBeReturned(mat.rows(), basis.rows());
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

  template <typename FieldType>
  Vector<FieldType> ColPivHouseholderQRSolver(Matrix<FieldType> A, 
                                              Vector<FieldType> b)
  { return A.colPivHouseholderQr().solve(b); }

  template <typename FieldType>
  Vector<FieldType> FullPivLUSolver(Matrix<FieldType> A, 
                                    Vector<FieldType> b)
  { return A.fullPivLu().solve(b); }

}

