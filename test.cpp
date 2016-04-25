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
	MatrixXf a(3, 2);
	a << 1, 2,
		 3, 4,
		 5, 6;
	MatrixXf b = a.rowwise() - a.colwise().mean();
	VectorXf c = stdev(a);
	cout << c << endl;
	b = b.array().rowwise() / c.transpose().array();
	cout << b << endl;
	//a.middleCols(0, 2).middleRows(1, 2) = b;
	//mat.conservativeResize(mat.rows() + 1, NoChange);
	//mat.row(mat.rows() - 1) = VectorXf::Ones(4);
	return 0;
}