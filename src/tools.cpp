#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	//declare and initialize rmse vector
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if((estimations.size() != ground_truth.size() )  || (estimations.size() == 0))
	{
		cout << "# of estimations does not correspond to # of ground truth samples " << endl;
		return rmse;
	}

	for(int i=0; i < estimations.size(); i++)
	{
		VectorXd res = estimations[i] - ground_truth[i];
		res = res.array()*res.array();
		rmse = rmse + res;
	}

	//mean value
	rmse = rmse/estimations.size();
	//squared root
	rmse = rmse.array().sqrt();

	return rmse;
}
