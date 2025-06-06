#pragma once

#include <Eigen/Dense>
// #include <mip/data_sensor.hpp>

class VelocityKalman3D {
public:
    VelocityKalman3D(double dtK = 0.01, double accel_var = 0.1, double meas_var = 0.5)
        : dtK(dtK), x(Eigen::Vector3d::Zero()), P(Eigen::Matrix3d::Identity())
    {
        A = Eigen::Matrix3d::Identity();
        B = Eigen::Matrix3d::Identity() * dtK;
        H = Eigen::Matrix3d::Identity();
        Q = Eigen::Matrix3d::Identity() * accel_var * dtK * dtK;
        R = Eigen::Matrix3d::Identity() * meas_var;
    }

    void step(
        double simcounter,
        const Eigen::Vector3d& accel_scaled,
        const Eigen::Matrix3d& R_body_to_world,//mip::data_sensor::CompOrientationMatrix& orientation_matrix,
        const double raw_velocity[3]
    ) {
        // Convert orientation matrix to Eigen
        // Eigen::Matrix3d R_body_to_world;
        // for (int i = 0; i < 3; ++i)
        //     for (int j = 0; j < 3; ++j)
        //         R_body_to_world(i, j) = orientation_matrix.matrix[i][j];

        // Rotate accel to world frame
        Eigen::Vector3d raw_velocity_estimate = Eigen::Vector3d::Zero(); // Placeholder for raw velocity estimate
        raw_velocity_estimate(0) = raw_velocity[0]; // Replace with actual raw velocity estimate
        raw_velocity_estimate(1) = raw_velocity[1]; // Replace with actual raw velocity estimate
        raw_velocity_estimate(2) = raw_velocity[2]; // Replace with actual raw velocity estimate
        Eigen::Vector3d a_world = R_body_to_world * accel_scaled;
        Eigen::Vector3d a_world2 = R_body_to_world.transpose() * accel_scaled;
        // std::cout << "accel_scaled" << "\t" << accel_scaled.transpose() << std::endl;
        // std::cout << "a_world" << "\t" << a_world.transpose() << std::endl;
        // Remove gravity (assuming Z-up)
        Eigen::Vector3d gravity(0, 0, 9.80665);
        Eigen::Vector3d a_corrected = a_world - gravity;

        // a_corrected(0) = a_corrected(0) > xddot_thresh ? xddot_thresh : (a_corrected(0) < -xddot_thresh ? -xddot_thresh : a_corrected(0)); 
        // a_corrected(1) = a_corrected(1) > yddot_thresh ? yddot_thresh : (a_corrected(1) < -yddot_thresh ? -yddot_thresh : a_corrected(1));
        // a_corrected(2) = a_corrected(2) > zddot_thresh ? zddot_thresh : (a_corrected(2) < -zddot_thresh ? -zddot_thresh : a_corrected(2)); 

        // Prediction
        x = A * x + B * a_corrected;
        P = A * P * A.transpose() + Q;

        // Kalman Gain
        Eigen::Matrix3d K = P * H.transpose() * (H * P * H.transpose() + R).inverse();

        // Update
        x = x + K * (raw_velocity_estimate - H * x);
        P = (Eigen::Matrix3d::Identity() - K * H) * P;

        // std::cout << simcounter << "\t" << a_corrected(0) << "\t" << a_corrected(1) << "\t" << a_corrected(2) << "\t" << raw_velocity[0] << "\t" << raw_velocity[1] << "\t" << raw_velocity[2] << "\t" 
        //                 << x(0) << "\t" << x(1) << "\t" << x(2) << std::endl;
        // std::cout << simcounter << "\t" << accel_scaled(0) << "\t" << accel_scaled(1) << "\t" << accel_scaled(2) << "\t" << a_world(0) << "\t" << a_world(1) << "\t" << a_world(2) << "\t" 
        //                 << a_world2(0) << "\t" << a_world2(1) << "\t" << a_world2(2) << std::endl;
    }

    Eigen::Vector3d getVelocity() const { return x; }

private:
    double dtK;
    Eigen::Vector3d x;
    Eigen::Matrix3d A, B, H, Q, R, P;
    double xddot_thresh = 10.0;
    double yddot_thresh = 10.0;
    double zddot_thresh = 10.0;
};

