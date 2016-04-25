#include "model.h"
#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

model::model() {}

model::model(const vector<int> &psize) {
	theta = MatrixXf::Zero(psize[0], psize[1]);
}

model::model(const vector<int> &psize, const vector<int> &tsize, const vector<int> &lsize) {
	theta = MatrixXf::Zero(psize[0], psize[1]);
	// Add ones for training data set
	x     = MatrixXf::Zero(tsize[0] + 1, tsize[1]);
	x.row(0).array().setConstant(1.0f);
	y     = MatrixXf::Zero(lsize[0], lsize[1]);
}

model::~model() {
	theta.resize(0, 0);
	x.resize(0, 0);
	y.resize(0, 0);
}

void model::getCost() {
	for (int i = 0; i < cost.size(); i++)
		cout << cost[i] << endl;
}

VectorXf stdev(const MatrixXf &mat) {
	return (((mat.rowwise() - mat.colwise().mean()).array().square() / mat.rows()).colwise().sum().array().sqrt()).row(0);
}

void model::initialize(const char *file, const char *mode) {
	/* mode = 'text' uses read function to parse text input file.
	   mode = 'binary' uses load function to load in the binary input file.
	*/
	MatrixXf mat;
	if (!strcmp(mode, "text"))
		mat = read(file);
	else if (!strcmp(mode, "binary"))
		mat = load(file);
	int rows = mat.rows(), cols = mat.cols();
	x.resize(rows, cols);
	x << MatrixXf::Ones(rows, 1), mat.leftCols(cols - 1);
	y.resize(rows, 1);
	y = mat.col(cols - 1);
	theta = MatrixXf::Zero(x.cols(), y.cols());
}

MatrixXf model::read(const char *file) {
	int cols = 0, rows = 0;
	vector<float> buff;

	ifstream infile;
	infile.open(file);

	while (!infile.eof()) {
		string line;
		getline(infile, line, '\n');
		int temp_cols = 0;
		float num;
		stringstream stream(line);
		while (stream >> num) {
			buff.push_back(num);
			temp_cols++;

			if (stream.peek() == ',')
				stream.ignore();
		}

		if (temp_cols == 0)
			continue;

		if (cols == 0)
			cols = temp_cols;

		rows++;
	}

	infile.close();

	MatrixXf mat(rows, cols);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			mat(i, j) = buff[cols * i + j];

	return mat;
}

MatrixXf model::load(const char *file) {
	ifstream infile(file, ios::in | ios::binary);
	typename MatrixXf::Index rows = 0, cols = 0;
	infile.read((char*) (&rows), sizeof(typename MatrixXf::Index));
	infile.read((char*) (&cols), sizeof(typename MatrixXf::Index));
	MatrixXf mat(rows, cols);
	infile.read((char*) mat.data(), rows * cols * sizeof(typename MatrixXf::Scalar));
	infile.close();
	return mat;
}

void model::save(const char *file, const MatrixXf &mat) {
	ofstream outfile(file, ios::out | ios::binary | ios::trunc);
	typename MatrixXf::Index rows = mat.rows(), cols = mat.cols();
	outfile.write((char*) (&rows), sizeof(typename MatrixXf::Index));
	outfile.write((char*) (&cols), sizeof(typename MatrixXf::Index));
	outfile.write((char*) mat.data(), rows * cols * sizeof(MatrixXf::Scalar));
	outfile.close();
}

void model::featureNormalize() {
	x = (x.rowwise() - x.colwise().mean()).array().rowwise() / (stdev(x)).transpose().array();
}

void model::LinearRegression(const float &alpha, const int &iter) {
	// Training set size
	int m = x.rows();
	// Feature size
	int n = x.cols();
	// Hypothesis - y
	MatrixXf temp;
	// Linear regression
	for (int i = 0; i < iter; i++) {
		temp = x * theta - y;
		// Cost function
		cost.push_back((temp.array().square()).sum() / (2 * m));
		// Update weigth
		theta.noalias() -= (alpha / m) * (x.adjoint() * temp);
	}
}

MatrixXf model::sigmoid(const MatrixXf &z) {
	return (1 / (1 + (-z).array().exp()));
}