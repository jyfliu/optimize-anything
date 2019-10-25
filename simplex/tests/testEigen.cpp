#include "../Eigen/Dense"
#include <iostream>

int main() {
    Eigen::MatrixXd mat(2, 2);
    mat.row(0) = Eigen::RowVector2d(6, 9);
    mat.row(1) = Eigen::RowVector2d(4, 2);
    std::cout << mat << std::endl;
}
