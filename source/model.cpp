#include "model.h"
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

/* Constructor */
model::model(const pair<int, int> &psize) {
	theta = MatrixXf::Zero(psize.first, psize.second);
}

model::model(const pair<int, int> &psize, const pair<int, int> &tsize, const pair<int, int> &lsize) {
	theta = MatrixXf::Zero(psize.first, psize.second);
	x     = MatrixXf::Zero(tsize.first, tsize.second);
	y     = MatrixXf::Zero(lsize.first, lsize.second);
}

model::model(const char *file, const char *mode) {
	/* 
	   mode = 'text' uses read function to parse text input file.
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
	x_mean.resize(1, cols);
	x_std.resize(cols);
}

/* Destructor */
model::~model() {
	theta.resize(0, 0);
	x.resize(0, 0);
	y.resize(0, 0);
}

/* Getter */
void model::getCost() {
	for (int i = 0; i < cost.size(); i++)
		cout << cost[i] << endl;
}

void model::getTheta() {
	cout << theta << endl;
}

/* Numerical computations outside API */
VectorXf stdev(const MatrixXf &mat) {
	return (((mat.rowwise() - mat.colwise().mean()).array().square() / mat.rows()).colwise().sum().array().sqrt()).row(0);
}

/* Initialization after construction */
void model::initialize(const char *file, const char *mode) {
	/* 
	   mode = 'text' uses read function to parse text input file.
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
	x_mean.resize(1, cols);
	x_std.resize(cols);
}

/* I/O */
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

/* Machine Learning */
void model::featureNormalize(MatrixXf &mat) {
	x_mean = mat.colwise().mean();
	x_std = stdev(mat);
	mat = (mat.rowwise() - x_mean.transpose()).array().rowwise() / x_std.transpose().array();
	mat.col(0).array().setConstant(1.0f);
}

void model::LinearRegression(const float &alpha, const int &iter, const float &lambda, const bool &useFeatureNormalize) {
	if (useFeatureNormalize)
		featureNormalize(x);
	gradientDescent(alpha, iter, lambda, "Linear");
}

void model::LogisticRegression(const float &alpha, const int &iter, const float &lambda, const bool &useFeatureNormalize) {
	if (useFeatureNormalize)
		featureNormalize(x);
	// Gradient Descent Solver
	//gradientDescent(alpha, iter, lambda, "Logistic");
	BFGS(iter, lambda);
}

MatrixXf model::sigmoid(const MatrixXf &z) {
	return (1 / (1 + (-z).array().exp()));
}

/* Solver */
MatrixXf model::gradientLinear(const float &alpha, const float &lambda) {
	// Training set size
	int m = x.rows();
	// Feature size
	int n = theta.rows();
	// Regularization term
	float reg = 0;
	// Hypothesis
	MatrixXf h = x * theta - y;
	// Gradient, same size as theta
	MatrixXf grad = x.adjoint() * h / m;
	// Cost function
	if (lambda != 0.0)	// Need to be verified
		reg = (2 * lambda / m) * (theta.bottomRows(n - 1).adjoint() * theta.bottomRows(n - 1)).array().sum();
	// Cost function
	cost.push_back(reg + (h.array().square()).sum() / (2 * m));
	if (lambda != 0.0)
		grad.bottomRows(n - 1).noalias() += lambda * theta.bottomRows(n - 1) / m;
	
	theta.noalias() -= alpha * grad;
	return grad;
}

MatrixXf model::gradientLogistic(const float &alpha, const float &lambda) {
	// Training set size
	int m = x.rows();
	// Feature size
	int n = theta.rows();
	// Regularization term
	float reg = 0;
	// Hypothesis
	MatrixXf h;
	// Gradient, same size as theta
	MatrixXf grad;
	// Data Holder
	MatrixXf H(2 * m, theta.cols()), Y(2 * m, y.cols());
	// Calculating cost
	h = sigmoid(x * theta);
	H << -h.array().log(),
		 -(1 - h.array()).array().log();

	h.noalias() -= y;

	Y << y,
		 1 - y.array();
	// Cost function
	if (lambda != 0.0)
		reg = (2 * lambda / m) * (theta.bottomRows(n - 1).adjoint() * theta.bottomRows(n - 1)).array().sum();
	
	cost.push_back(reg + (H.adjoint() * Y / m).array().sum());
	
	grad = x.adjoint() * h / m;

	if (lambda != 0.0)
		grad.bottomRows(n - 1).noalias() += lambda * theta.bottomRows(n - 1) / m;

	theta.noalias() -= alpha * grad;
	return grad;
}

void model::gradientDescent(const float &alpha, const int &iter, const float &lambda, const char *mode) {
	// Regression loop
	for (int i = 0; i < iter; i++) {
		if (!strcmp(mode, "Logistic"))
			(void)gradientLogistic(alpha, lambda);

		else if (!strcmp(mode, "Linear"))
			(void)gradientLinear(alpha, lambda);
	}
}

void model::BFGS(const int &maxiter, const float &lambda) {
	// Feature size
	int n = theta.rows();
	// Alpha
	float alpha = 15;
	// Output and Output update
	MatrixXf Fc, Fc_new, df, d;
	// Hessian Matrix
	MatrixXf Hessian = MatrixXf::Identity(n, n);
	// Step difference
	MatrixXf s, q;
	// Initial Step
	// Define Gradient
	MatrixXf grad, grad_new;
	// Gradient
	grad = gradientLogistic(alpha, lambda);
	// Start of BFGS
	s = Hessian.ldlt().solve(-alpha * grad);	// theta_n+1 - theta_n
	// Update theta
	// Get current function output
	Fc = sigmoid(x * theta);
	d  = -Hessian * grad;
	df = grad.adjoint() * d;
	Fc_new = sigmoid(x * (theta + d));
	// Line search
	while ((Fc_new.array() > (Fc + alpha * df).array()).any()) {
		alpha *= 0.7; 
		Fc_new = sigmoid(x * (theta + alpha * d));
		//cout << "alpha: " << alpha << endl;
	}
	// Loop
	for (int i = 1; i < maxiter; i++) {
		cout << "alpha: " << alpha << endl;
		// Gradient
		grad_new = gradientLogistic(alpha, lambda);
		if (cost.rbegin()[1] - cost.rbegin()[0] < 1e-8) {
			cout << "Early stop criteria met!" << endl;
			cout << "Loop time: " << i << endl;
			cout << (cost.rbegin()[1] - cost.rbegin()[0]) << endl;
			break;
		}
		// Gradient difference between step
		q = grad_new - grad;
		//cout << "Grad: " << endl << grad << endl;
		//cout << "Grad_New: " << endl << grad_new << endl;
		grad = grad_new;
		// Calculate Hessian Matrix
		//cout << "Strictly Positive? " << q.adjoint() * s << endl;
		Hessian.noalias() += (q * q.adjoint()) / (q.adjoint() * s)(0) - (Hessian * s * s.adjoint() * Hessian) / (s.adjoint() * Hessian * s)(0);
		s = Hessian.ldlt().solve(-alpha * grad);	// theta_n+1 - theta_n
		// Line search
		alpha = 15;
		Fc = sigmoid(x * theta);
		d  = -Hessian * grad;
		df = grad.adjoint() * d;
		Fc_new = sigmoid(x * (theta + d));
		while ((Fc_new.array() > (Fc + alpha * df).array()).any()) {
			alpha *= 0.7; 
			Fc_new = sigmoid(x * (theta + alpha * d));
		}

		// End of line search
	}
}

/* Output */
void model::predict(MatrixXf &mat, const bool &useFeatureNormalize) {
	if (useFeatureNormalize) {
		mat = (mat.rowwise() - x_mean.transpose()).array().rowwise() / x_std.transpose().array();
		mat.col(0).array().setConstant(1.0f);
	}
	cout << sigmoid(mat * theta) << endl;
}

/* Timer */
void model::timer(const char *mode) {
	/*
	   mode = start : Setting the time_start to now
	   mode = stop  : Printing the time duration since last start
	*/
	if (!strcmp(mode, "start"))
		time_start = chrono::high_resolution_clock::now();
	else if (!strcmp(mode, "stop"))
		cout << "Elapsed time is " \
	  		 << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - time_start).count() \
	  		 << " ms." << endl;
	else
		throw invalid_argument("Bad mode!");
}