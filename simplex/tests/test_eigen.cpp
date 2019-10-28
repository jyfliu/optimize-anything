#include "../Eigen/Dense"
#include <iostream>

void t1() {
  Eigen::MatrixXd mat(2, 2);
  mat.row(0) = Eigen::RowVector2d(6, 9);
  mat.row(1) = Eigen::RowVector2d(4, 2);
  //std::cout << mat << std::endl;
}

void t2() {
  Eigen::MatrixXd mat(3, 3);
  mat << 3, 1, 3, 2, 2, 1, -2*3+5*2, -2*1+5*2, -2*3+5*1;
  std::cout << Eigen::MatrixXd(mat.colPivHouseholderQr().matrixQ()) << std::endl;
  std::cout << Eigen::MatrixXd(mat.colPivHouseholderQr().matrixQR().triangularView<Eigen::Upper>()) << std::endl;
}

int main() {
  t1();
  t2();

  return 0;
}
