#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  // trigger initialization
  is_initialized_ = false;

  // dt
  time_us_ = 0.0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  ///* State dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  ///* Augmented state dimension
  n_aug_ = 7;

  //Sigma points scaling rel. to mean
  lambda_ = 3 - n_aug_;

  //Matrix including predicted Sigma Points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // weights
  weights_ = VectorXd(2*n_aug_+1);

  // radar specific measurement noise covariance matrix
  Noise_radar_ = MatrixXd(3,3);
  Noise_radar_ << std_radr_*std_radr_, 0, 0,
  	  		 	  0, std_radphi_*std_radphi_, 0,
  				  0, 0,std_radrd_*std_radrd_;

  // lidar specific measurement noise covariance matrix
  Noise_lidar_ = MatrixXd(2,2);
  Noise_lidar_	 << std_laspx_* std_laspx_, 0,
		  	  	  	  0, std_laspy_* std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

/*****************************************************************************
   *  Initialization
   ****************************************************************************/
	if (!is_initialized_)
	{

	cout << "UKF: INIT " << endl;

	/** Set initial values to Identity */
	P_ 	<< 	  1, 0, 0, 0, 0,
	          0, 1, 0, 0, 0,
	          0, 0, 1, 0, 0,
	          0, 0, 0, 1, 0,
	          0, 0, 0, 0, 1;

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
	  /**
	  Convert radar from polar to cartesian coordinates
	  */
			float rho = meas_package.raw_measurements_[0];
			float phi = meas_package.raw_measurements_[1];

			float px = rho * cos(phi);
			float py = rho * sin(phi);

			// update state vector with first values.
			// do not use vx for now (treat as hidden state)
			x_ << px, py, 0, 0, 0;
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
			/** this is simple as no conversion necessary. Set the velocities to 0
				as these represent the hidden part (that is to be observed) **/
				x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
	}

	// store first timestamp
	time_us_ = meas_package.timestamp_;


	/** Set weights */
	weights_ = SetWeights(lambda_, n_aug_);

	// done initializing, no need to predict or update
	is_initialized_ = true;
	cout << "UKF: INIT DONE " << endl;
	return;
	}

/*****************************************************************************
   *  Normal Run
   ****************************************************************************/
	float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;

	Predict(dt);
	Update(meas_package);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Predict(double delta_t) {

	/** Generate augmented X matrix */
	MatrixXd Xsig_aug = GenerateAugmentedSigmaPoints(n_aug_, lambda_, std_a_, std_yawdd_, x_, P_);

	/** Predict Sigma Points */
	Xsig_pred_ = PredictSigmaPoints(n_x_, n_aug_, Xsig_aug, delta_t);

	/** Predict Mean */
	x_ = PredictMean(n_aug_, weights_, Xsig_pred_);

	/** Predict Covariance */
	P_ = PredictCovariance(n_aug_,weights_, Xsig_pred_, x_);
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::Update(MeasurementPackage meas_package) {

	int n_z;
	VectorXd z;
	MatrixXd Zsig;
	MatrixXd Noise;
	bool normalizeAngles;

	if ( (meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_ == true))
	{/* Radar updates */

		// radar specific dimensions
		n_z = 3;
		// radar noise
		Noise = Noise_radar_;
		normalizeAngles = true;

		Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
		//transform sigma points into measurement space
		for (int i = 0; i < 2 * n_aug_ + 1; i++)
		{
			// extract values
			double p_x 	= Xsig_pred_(0,i);
			double p_y 	= Xsig_pred_(1,i);
			double v  	= Xsig_pred_(2,i);
			double yaw 	= Xsig_pred_(3,i);

			double v1 	= cos(yaw)*v;
			double v2 	= sin(yaw)*v;
			// measurement model
			Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
			Zsig(1,i) = atan2(p_y,p_x);                                 //phi
			Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
		}
	}
	else if( (meas_package.sensor_type_ == MeasurementPackage::LASER ) && (use_laser_ == true) )
	{/* Lidar updates */

		// lidar specific dimensions
		n_z = 2;
		// lidar noise
		Noise = Noise_lidar_;
		normalizeAngles = false;

		Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
		//transform sigma points into measurement space
		for (int i = 0; i < 2 * n_aug_ + 1; i++)
		{
			// extract values
			double p_x = Xsig_pred_(0,i);
			double p_y = Xsig_pred_(1,i);
			// measurement model
			Zsig(0,i) = p_x;
			Zsig(1,i) = p_y;
		}
	}
	else
	{/* Ground truth values. Do not evaluate. */
		return;
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++)
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	MatrixXd S = CalculateInnovationCovarianceMatrix(n_z, n_aug_, weights_, Zsig, z_pred, Noise, normalizeAngles);

	MatrixXd K = CalculateKalmanGain(n_x_, n_z, n_aug_, Zsig, Xsig_pred_, weights_, S, z_pred, x_, normalizeAngles);

	z = VectorXd(n_z);
	z =  meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred;

	if(normalizeAngles)
	{
		//angle normalization
		while (z_diff(1)> M_PI){z_diff(1)-=2.*M_PI;}
		while (z_diff(1)<-M_PI){z_diff(1)+=2.*M_PI;}
	}

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
}



/*************************************************************
** HELPER FUNCTIONS FOR PREDICTION STEP						*
**************************************************************/
/**
 * Generate augmented Sigma Points
 * @param {int} n_aug the dimension of the augmented matrix
 * @param {double} lambda the distance of the Sigma Points relative to the mean
 * @param {double} std_a the process noise standard deviation longitudinal acceleration
 * @param {double} std_yawdd the process noise standard deviation yaw acceleration
 * @param {VectorXd} x the state vector
 * @param {MatrixXd} P the covariance matrix */
MatrixXd UKF::GenerateAugmentedSigmaPoints(int n_aug, double lambda, double std_a, double std_yawdd, VectorXd x, MatrixXd P){

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug);
	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug, n_aug);
	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
	//create augmented mean state
	x_aug.head(5) = x;
	x_aug(5) = 0;
	x_aug(6) = 0;
	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P;
	P_aug(5,5) = std_a*std_a;
	P_aug(6,6) = std_yawdd*std_yawdd;
	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();
	//create augmented sigma points
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i< n_aug; i++)
	{
		Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
		Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
	}
	return Xsig_aug;
}

/** Prediction of Sigma Points
 * @param {int} n_x the dimension of the state vector
 * @param {int} n_aug the dimension of the augmented matrix
 * @param {MatrixXd} Xsig_aug the augmented Sigma Points
 * @param {double} delta_t the time step */
MatrixXd UKF::PredictSigmaPoints(int n_x, int n_aug, MatrixXd Xsig_aug, double delta_t){

	MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

	//predict sigma points
	for (int i = 0; i< 2*n_aug+1; i++)
	{
		//extract values for better readability
		double p_x 		= Xsig_aug(0,i);
		double p_y 		= Xsig_aug(1,i);
		double v 		= Xsig_aug(2,i);
		double yaw 		= Xsig_aug(3,i);
		double yawd 		= Xsig_aug(4,i);
		double nu_a 		= Xsig_aug(5,i);
		double nu_yawdd 	= Xsig_aug(6,i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred(0,i) = px_p;
		Xsig_pred(1,i) = py_p;
		Xsig_pred(2,i) = v_p;
		Xsig_pred(3,i) = yaw_p;
		Xsig_pred(4,i) = yawd_p;
	}
	return Xsig_pred;
}

/**Set weights
 * @param {double} lambda the distance of the Sigma Points relative to the mean
 * @param {int} n_aug the dimension of the augmented matrix */
VectorXd UKF::SetWeights(double lambda, int n_aug){

	VectorXd weights = VectorXd(2*n_aug+1);
	double weight_0 = lambda/(lambda+n_aug);
	weights(0) = weight_0;
	for (int i=1; i<2*n_aug+1; i++)
	{
		double weight = 0.5/(n_aug+lambda);
		weights(i) = weight;
	}
	return weights;
}

/** Predict state mean
 * @param {int} n_aug the dimension of the augmented matrix
 * @param {MatrixXd} Xsig_pred the predicted sigma points
 * @param {VectorXd} weights the weights sigma points based on the lambda value  */
MatrixXd UKF::PredictMean(int n_aug, VectorXd weights, MatrixXd Xsig_pred){

	VectorXd x = VectorXd(5);
	 x.fill(0.0);
	 for (int i = 0; i < 2 * n_aug + 1; i++)
	 {  //iterate over sigma points
		 x = x+ weights(i) * Xsig_pred.col(i);
	 }
	return x;
}

/** Predict covariance
 * @param {int} n_aug the dimension of the augmented matrix
 * @param {VectorXd} weights the weights sigma points based on the lambda value
 * @param {MatrixXd} Xsig_pred the predicted sigma points
 * @param {VectorXd} x the state */
MatrixXd UKF::PredictCovariance(int n_aug, VectorXd weights, MatrixXd Xsig_pred, VectorXd x){

	//predicted state covariance matrix
	MatrixXd P = MatrixXd(5, 5);
	//predicted state covariance matrix
	P.fill(0.0);
	for (int i = 0; i < 2 * n_aug + 1; i++)
	{  //iterate over sigma points
	    // state difference
		VectorXd x_diff = Xsig_pred.col(i) - x;
	    //angle normalization
	    while (x_diff(3)> M_PI){x_diff(3)-=2.*M_PI;}
	    while (x_diff(3)< -M_PI){x_diff(3)+=2.*M_PI;}
	    P = P + weights(i) * x_diff * x_diff.transpose() ;
	  }
	return P;
}


/*************************************************************
** HELPER FUNCTIONS FOR UPDATE STEP  						*
**************************************************************/

/* Calculation of the innovation covariance Matrix S
 * @param {int} n_z the dimension of the measurement vector
 * @param {int} n_aug the dimension of the augmented matrix
 * @param {VectorXd} weights the weights sigma points based on the lambda value
 * @param {MatrixXd} Zsig the measurement matrix
 * @param {VectorXd} z_pred the predicted state
 * @param {MatrixXd} R the noise covariance matrix
 * @param {bool} normalizeAngles */
MatrixXd UKF::CalculateInnovationCovarianceMatrix(int n_z, int n_aug, VectorXd weights, MatrixXd Zsig, VectorXd z_pred, MatrixXd Noise, bool normalizeAngles){

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug + 1; i++)
	{  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		if (normalizeAngles)
		{
			//angle normalization
			while (z_diff(1)> M_PI){z_diff(1)-=2.*M_PI;}
			while (z_diff(1)<-M_PI){z_diff(1)+=2.*M_PI;}
		}
		S = S + weights(i) * z_diff * z_diff.transpose();
	}
	S = S + Noise;
	return S;
}

/**
 * Calculation of the Kalman Gain
 * @param {int} n_x the dimension of the state vector
 * @param {int} n_aug the dimension of the augmented matrix
 * @param {MatrixXd} Zsig the measurement matrix
 * @param {MatrixXd} Xsig_pred the predicted Sigma Points
 * @param {VectorXd} weights the weights sigma points based on the lambda value
 * @param {MatrixXd} S the innovation covariance matrix
 * @param {VectorXd} z_pred the predicted state
 * @param {VectorXd} x the state vector
 * @param {bool} normalizeAngles */
MatrixXd UKF::CalculateKalmanGain(int n_x, int n_z, int n_aug, MatrixXd Zsig, MatrixXd Xsig_pred, VectorXd weights, MatrixXd S, VectorXd z_pred, VectorXd x, bool normalizeAngles){

	MatrixXd Tc = MatrixXd(n_x, n_z);
	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug + 1; i++)
	{
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		if(normalizeAngles)
		{
			//angle normalization
			while (z_diff(1)> M_PI){ z_diff(1)-=2.*M_PI;}
			while (z_diff(1)< -M_PI){ z_diff(1)+=2.*M_PI;}
		}
		// state difference
		VectorXd x_diff = Xsig_pred.col(i) - x;

		if(normalizeAngles)
		{
			//angle normalization
			while (x_diff(3)> M_PI) {x_diff(3)-=2.*M_PI;} // issue here
			while (x_diff(3)<-M_PI) {x_diff(3)+=2.*M_PI;}
		}

		Tc = Tc + weights(i) * x_diff * z_diff.transpose();
	}
	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	return K;
}
