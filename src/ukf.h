#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:


  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  // radar specific measurement noise covariance matrix
  MatrixXd Noise_radar_;
  // radar specific measurement noise covariance matrix
  MatrixXd Noise_lidar_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Predict(double delta_t);

  /**
     * Updates the state and the state covariance matrix
     * @param meas_package The measurement at k+1
     */
  void Update(MeasurementPackage meas_package);

  /**
   * Generate augmented Sigma Points
   * @param {int} n_aug the dimension of the augmented matrix
   * @param {int} lambda the distance of the Sigma Points relative to the mean
   * @param {int} std_a the process noise standard deviation longitudinal acceleration
   * @param {int} std_yawdd the process noise standard deviation yaw acceleration
   * @param {VectorXd} x the state vector
   * @param {MatrixXd} P the covariance matrix
   */
  MatrixXd GenerateAugmentedSigmaPoints(int n_aug, int lambda, int std_a, int std_yawdd, VectorXd x, MatrixXd P);

  /** Prediction of Sigma Points
   * @param {int} n_x the dimension of the state vector
   * @param {int} n_aug the dimension of the augmented matrix
   * @param {MatrixXd} Xsig_aug the augmented Sigma Points
   * @param {int} delta_t the time step
   */
  MatrixXd PredictSigmaPoints(int n_x, int n_aug, MatrixXd Xsig_aug, double delta_t);

  /**Set weights
   * @param {int} lambda the distance of the Sigma Points relative to the mean
   * @param {int} n_aug the dimension of the augmented matrix
   */
  VectorXd SetWeights(int lambda, int n_aug);

  /** Predict state mean
   * @param {int} n_aug the dimension of the augmented matrix
   * @param {MatrixXd} Xsig_pred the predicted sigma points
   * @param {VectorXd} weights the weights sigma points based on the lambda value
   */
  MatrixXd PredictMean(int n_aug, VectorXd weights, MatrixXd Xsig_pred);

  /** Predict covariance
   * @param {int} n_aug the dimension of the augmented matrix
   * @param {VectorXd} weights the weights sigma points based on the lambda value
   * @param {MatrixXd} Xsig_pred the predicted sigma points
   * @param {VectorXd} x the state
   */
  MatrixXd PredictCovariance(int n_aug, VectorXd weights, MatrixXd Xsig_pred, VectorXd x);

  /* Calculation of the innovation covariance Matrix S
   * @param {int} n_z the dimension of the measurement vector
   * @param {int} n_aug the dimension of the augmented matrix
   * @param {VectorXd} weights the weights sigma points based on the lambda value
   * @param {MatrixXd} Zsig the measurement matrix
   * @param {VectorXd} z_pred the predicted state
   * @param {MatrixXd} R the noise covariance matrix
   */
  MatrixXd CalculateInnovationCovarianceMatrix(int n_z, int n_aug, VectorXd weights, MatrixXd Zsig, VectorXd z_pred, MatrixXd Noise, bool normalizeAngles);

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
   */
  MatrixXd CalculateKalmanGain(int n_x, int n_z, int n_aug, MatrixXd Zsig, MatrixXd Xsig_pred, VectorXd weights, MatrixXd S, VectorXd z_pred, VectorXd x, bool normalizeAngles);

};

#endif /* UKF_H */
