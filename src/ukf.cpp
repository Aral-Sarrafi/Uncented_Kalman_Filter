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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.
  
  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  previous_timestamp = 0;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if (!is_initialized_) {
		// Initiate the state vector
		x_ << 0, 0, 0, 0, 0;

		// Initiate the Covariance Matrix
		P_.fill(0.0);
		for (int i = 0; i < x_.size(); i++)
		{
			P_(i, i) = 1;
		}

		// Use the first measurement info as the first state vector
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			x_(0) = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
			x_(1) = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

			x_(0) = meas_package.raw_measurements_[0];
			x_(0) = meas_package.raw_measurements_[1];
		}

		previous_timestamp = meas_package.timestamp_;
		is_initialized_ = true;

		return;
	}
	else
	{
		double delta_t = (meas_package.timestamp_ - previous_timestamp) / 1000000.0;
		previous_timestamp = meas_package.timestamp_;

		Prediction(delta_t);

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			UpdateRadar(meas_package);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

			UpdateLidar(meas_package);

		}

	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
	// Generate the augmented state vector
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;

	// Generate the augmented covariance matrix
	MatrixXd Q = MatrixXd(2, 2);
	Q.fill(0.0);
	Q(0, 0) = std_a_ * std_a_;
	Q(1, 1) = std_yawdd_ * std_yawdd_;

	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug.bottomRightCorner(2, 2) = Q;

	// Square root of P_aug matrix
	MatrixXd A = P_aug.llt().matrixL();

	// Generate the augmented sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		if (i == 0)
		{
			Xsig_aug.col(i) = x_aug;
		}
		else if (i <= n_aug_ && i > 0)
		{
			Xsig_aug.col(i) = x_aug + sqrt(lambda_ + n_aug_)*A.col(i - 1);

		}
		else if (i > n_aug_) {
			Xsig_aug.col(i) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i - n_aug_ - 1);
		}

	}

	// Sigma point prediction
	double delta_t2 = delta_t * delta_t;



	for (int i = 0; i< 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		}
		else {
			px_p = p_x + v * delta_t*cos(yaw);
			py_p = p_y + v * delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd * delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

	// Predic the mean and the covariance matrix
	//create vector for predicted state
	VectorXd x = VectorXd(n_x_);

	//create covariance matrix for prediction
	MatrixXd P = MatrixXd(n_x_, n_x_);


	for (int i = 1; i < 2 * n_aug_ + 1; i++)
	{
		weights_(i) = 0.5 / (lambda_ + n_aug_);
	}
	//predict state mean
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x += weights_(i)*Xsig_pred_.col(i);
	}

	//predict state covariance matrix
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		P += weights_(i)*(Xsig_pred_.col(i) - x)*(Xsig_pred_.col(i) - x).transpose();
	}
	
	// Update the state vectors and covariance matrix
	x_ = x;
	P_ = P;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	float px = meas_package.raw_measurements_[0];
	float py = meas_package.raw_measurements_[1];


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	float rho = meas_package.raw_measurements_[0];
	float phi = meas_package.raw_measurements_[1];
	float rho_dot = meas_package.raw_measurements_[2];



}
