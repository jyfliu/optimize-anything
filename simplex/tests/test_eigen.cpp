#include "../Eigen/Dense"
#include <iostream>

// sandbox
void t1() {
  Eigen::MatrixXd mat(2, 2);
  mat.row(0) = Eigen::RowVector2d(6, 9);
  mat.row(1) = Eigen::RowVector2d(4, 2);
  std::cout << mat << std::endl;
}

void t2() {
  Eigen::MatrixXd mat(3, 3);
  mat << 3, 1, 3, 2, 2, 1, -2*3+5*2, -2*1+5*2, -2*3+5*1;
  std::cout << Eigen::MatrixXd(mat.colPivHouseholderQr().matrixQ()) << std::endl;
  std::cout << Eigen::MatrixXd(mat.colPivHouseholderQr().matrixQR().triangularView<Eigen::Upper>()) << std::endl;
}

void t3() {
  using std::cout;
  using std::endl;
  Eigen::MatrixXd m(3, 3);
  m << 3, 1, 3, 2, 2, 1, -2*3+5*2, -2*1+5*2, -2*3+5*1;
  cout << "Here is the matrix m:" << endl << m << endl;
  auto lu = m.fullPivLu();
  cout << "Here is, up to permutations, its LU decomposition matrix:"
       << endl << lu.matrixLU() << endl;
  cout << "Here is the L part:" << endl;
  Eigen::MatrixXd l = Eigen::MatrixXd::Identity(3, 3);
  l.triangularView<Eigen::StrictlyLower>() = lu.matrixLU();
  cout << l << endl;
  cout << "Here is the U part:" << endl;
  Eigen::MatrixXd u = lu.matrixLU().triangularView<Eigen::Upper>();
  cout << u << endl;
  cout << "Let us now reconstruct the original matrix m:" << endl;
  cout << lu.permutationP().inverse() * l * u * lu.permutationQ().inverse() << endl;
}

int main() {
  t3();
  return 0;
}
