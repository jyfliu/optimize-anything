#include "../../src/utils/eigen_utils.hpp"

#include <iostream>

void test_subRows() {
#define ASSERT(x) assert(x and "test_subRows ")
  Eigen::MatrixXd mat(5, 5);
  mat << 1, 2, 3, 4, 5,
         6, 7, 8, 9, 10,
         11, 12, 13, 14, 15,
         16, 17, 18, 19, 20,
         21, 22, 23, 24, 25;
  Eigen::VectorXi vec(3);
  vec << 0, 2, 3;
  Eigen::MatrixXd res(3, 5);
  res << 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
  ASSERT(simplex::subRows(mat, vec) == res);
#undef ASSERT
}

int main() {
  test_subRows();
  std::cout << "All tests passed." << std::endl;
  return 0;
}
