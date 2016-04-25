#include "Eigen/Dense"
#include <vector>

using namespace std;
using namespace Eigen;

class model {
public:
	/* Constructor and Destructor */
	model();
	model(const vector<int> &psize);
	model(const vector<int> &psize, const vector<int> &tsize, const vector<int> &lsize);
	~model();

	void getCost();

	void initialize(const char *file, const char *mode);

	/* I/O */
	MatrixXf read(const char *file);
	MatrixXf load(const char *file);
	void featureNormalize();
	void save(const char *file, const MatrixXf &mat);

	/* Regresssion */
	void LinearRegression(const float &alpha, const int &iter);
	MatrixXf sigmoid(const MatrixXf &z);

private:
	MatrixXf theta;
	MatrixXf x;
	MatrixXf y;
	vector<float> cost;
};