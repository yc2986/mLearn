#include <iostream>
#include <chrono>
#include "Eigen/Dense"
#include "model.h"

using namespace std;
using namespace Eigen;

int main() {
	float alpha;
	int iter;
	cout << "Please input the learning rate and iteration time!" << endl;
	cin >> alpha >> iter;
	vector<int> psize = {1, 1};
	model foo(psize);
	MatrixXf temp = foo.read("model/ex1data2.txt");
	foo.save("model/out2.csv", temp);
	foo.initialize("model/out2.csv", "binary");
	auto start = chrono::high_resolution_clock::now();
	foo.featureNormalize();
	foo.LinearRegression(alpha, iter);
	foo.getCost();
	auto dur = chrono::high_resolution_clock::now() - start;
	auto ms = chrono::duration_cast<chrono::milliseconds>(dur).count();
	cout << ms << "ms" << endl;
	return 0;
}