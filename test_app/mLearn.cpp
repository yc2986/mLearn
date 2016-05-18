#include <iostream>
#include <utility>
#include "Eigen/Dense"
#include "model.h"

using namespace std;
using namespace Eigen;

int main() {
	float alpha;
	int iter;
	MatrixXf score(1, 3);
	score << 1, 45, 85;
	cout << "Please input the learning rate and iteration time!" << endl;
	cin >> alpha >> iter;
	model foo("model/ex2data1.txt", "text");
	foo.timer("start");
	foo.LogisticRegression(alpha, iter, 0, true);
	foo.getCost();
	foo.predict(score, true);
	foo.timer("stop");
	return 0;
}