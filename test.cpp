#include <iostream>
#include <string>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;
/*
void push_back(MatrixXf &mat, VectorXf &vec) {
	int row = mat.rows();
	mat.conservativeResize(row + 1, Eigen::NoChange);
	mat.row(row) = vec;
}
*/

VectorXf stdev(MatrixXf mat) {
	return (((mat.rowwise() - mat.colwise().mean()).array().square() / mat.rows()).colwise().sum().array().sqrt()).row(0);
}

int main() {
	MatrixXf a(4, 2);
	VectorXf b(2);
	a << 1, 5,
		 2, 6,
		 3, 7,
		 4, 8;
	b = a.colwise().mean();
	a.rowwise() -= b.transpose();
	cout << a << endl;
	return 0;
}