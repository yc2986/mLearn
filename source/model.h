#include "Eigen/Dense"
#include <utility>
#include <vector>
#include <chrono>

using namespace std;
using namespace Eigen;

class model {
public:
	/* Constructor */
	model() {};
	model(const pair<int, int> &psize);
	model(const pair<int, int> &psize, const pair<int, int> &tsize, const pair<int, int> &lsize);
	model(const char *file, const char *mode);

	/* Destructor */
	~model();

	/* Getter */
	void getCost();
	void getTheta();

	/* Initialization after construction */
	void initialize(const char *file, const char *mode);

	/* I/O */
	MatrixXf read(const char *file);
	MatrixXf load(const char *file);
	void save(const char *file, const MatrixXf &mat);

	/* Regresssion */
	void LinearRegression(const float &alpha, const int &iter, const float &lambda = 0, const bool &useFeatureNormalize = false);
	void LogisticRegression(const float &alpha, const int &iter, const float &lambda = 0,const bool &useFeatureNormalize = false);
	MatrixXf sigmoid(const MatrixXf &z);

	/* Output */
	void predict(MatrixXf &mat, const bool &useFeatureNormalize = false);

	/* Timer */
	void timer(const char *mode);

private:
	MatrixXf theta;
	MatrixXf x;
	MatrixXf y;
	vector<float> cost;
	VectorXf x_mean;
	VectorXf x_std;
	chrono::high_resolution_clock::time_point time_start;

	void featureNormalize(MatrixXf &mat);
};