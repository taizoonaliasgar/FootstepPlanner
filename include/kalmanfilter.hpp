// vel_bias_kf.hpp
#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <limits>

template<int D = 3>
class VelBiasKF {
public:
    static_assert(D >= 1, "D must be >= 1");
    using VecD  = Eigen::Matrix<double, D, 1>;
    using MatD  = Eigen::Matrix<double, D, D>;
    using Vec2D = Eigen::Matrix<double, 2*D, 1>;
    using Mat2D = Eigen::Matrix<double, 2*D, 2*D>;

    VelBiasKF() { reset(); }

    // Reset to sane defaults
    void reset() {
        x_.setZero();
        P_.setZero();
        P_.template block<D,D>(0,0).setIdentity();  P_.template block<D,D>(0,0) *= 1e-1; // vel
        P_.template block<D,D>(D,D).setIdentity();  P_.template block<D,D>(D,D) *= 1e-2; // bias
        R_.setIdentity(); R_ *= 1e-3;
        scale_qv_ = 1.0;
        scale_qb_ = 1e-4;
        sigma_a2_ = 1e-4; // accel variance estimate (m/s^2)^2
    }

    // Initialize state
    void setInitial(const VecD& v0, const VecD& b0 = VecD::Zero()) {
        x_.template segment<D>(0) = v0;
        x_.template segment<D>(D) = b0;
    }

    // Tuning knobs
    void setMeasurementCovariance(const MatD& R) { R_ = R; }     // R ~ var(v_meas)
    void setProcessScales(double qv_scale, double qb_scale) { scale_qv_ = qv_scale; scale_qb_ = qb_scale; }
    void setAccelVariance(double sigma_a2) { sigma_a2_ = std::max(sigma_a2, 1e-12); } // avg accel var
    void settimestep(double dtKF) { dt_k = std::max(dtKF, 1e-6); }
    // One filtering step.
    // a_meas: gravity-compensated accel in world frame [m/s^2]
    // v_meas: velocity measurement (same frame) [m/s]
    // dt    : seconds (must be > 0)
    void step(VecD& a_meas, const VecD& v_meas, VecD* v_hat_out = nullptr, VecD* b_hat_out = nullptr) {
        // const double dt_k = (dtKF > 0.0 && std::isfinite(dtKF)) ? dtKF : std::numeric_limits<double>::epsilon();


        a_meas(0) = a_meas(0) > xddot_thresh ? xddot_thresh : (a_meas(0) < -xddot_thresh ? -xddot_thresh : a_meas(0));
        a_meas(1) = a_meas(1) > yddot_thresh ? yddot_thresh : (a_meas(1) < -yddot_thresh ? -yddot_thresh : a_meas(1));
        a_meas(2) = a_meas(2) > zddot_thresh ? zddot_thresh : (a_meas(2) < -zddot_thresh ? -zddot_thresh : a_meas(2));
        // Build F, G
        Mat2D F = Mat2D::Identity();
        F.template block<D,D>(0,D) = -MatD::Identity() * dt_k; // dv = ... - b*dt
        Eigen::Matrix<double, 2*D, D> G; G.setZero();
        G.template block<D,D>(0,0) = MatD::Identity() * dt_k;  // input = accel

        // Process noise Q (velocity jerk + bias random walk)
        const double qv = std::max(scale_qv_ * (sigma_a2_ * dt_k * dt_k), QV_FLOOR);
        const double qb = std::max(scale_qb_ * dt_k, QB_FLOOR);
        Mat2D Q = Mat2D::Zero();
        Q.template block<D,D>(0,0).setIdentity();  Q.template block<D,D>(0,0) *= qv;
        Q.template block<D,D>(D,D).setIdentity();  Q.template block<D,D>(D,D) *= qb;

        // Predict
        VecD u = a_meas;
        for (int i = 0; i < D; ++i) if (!std::isfinite(u(i))) u(i) = 0.0;
        Vec2D x_pred = F * x_ + G * u;
        Mat2D P_pred = F * P_ * F.transpose() + Q;

        // Update only if v_meas is finite
        bool meas_ok = true;
        for (int i = 0; i < D; ++i) if (!std::isfinite(v_meas(i))) { meas_ok = false; break; }

        if (meas_ok) {
            // H = [I_D  0]
            Eigen::Matrix<double, D, 2*D> H; H.setZero(); H.template block<D,D>(0,0).setIdentity();

            MatD S = (H * P_pred * H.transpose()) + R_;
            MatD Ssym = (S + S.transpose()) * 0.5;   // enforce symmetry

            // Robust solve for K = P_pred*H' * S^{-1}
            Eigen::LDLT<MatD> ldlt(Ssym);
            int tries = 0;
            while (ldlt.info() != Eigen::Success && tries < 5) {
                Ssym.diagonal().array() += 1e-9;     // jitter
                ldlt.compute(Ssym);
                ++tries;
            }

            Eigen::Matrix<double, 2*D, D> K;
            if (ldlt.info() == Eigen::Success) {
                const auto PHt = P_pred * H.transpose();
                // solve S * X = (PHt)^T  => X^T is K
                K = ldlt.solve(PHt.transpose()).transpose();
            } else {
                // last resort (should be rare)
                K = P_pred * H.transpose() * Ssym.inverse();
            }

            const VecD z = v_meas;
            const VecD innov = z - x_pred.template segment<D>(0);   // z - H*x_pred
            x_ = x_pred + K * innov;

            const Mat2D I = Mat2D::Identity();
            P_ = (I - K*H) * P_pred * (I - K*H).transpose() + K * R_ * K.transpose(); // Joseph form
        } else {
            // predict-only step
            x_ = x_pred;
            P_ = P_pred + Mat2D::Identity() * 1e-12;
        }

        // Symmetrize covariance
        P_ = (P_ + P_.transpose()) * 0.5;

        if (v_hat_out) *v_hat_out = x_.template segment<D>(0);
        if (b_hat_out) *b_hat_out = x_.template segment<D>(D);
    }

    // Accessors
    VecD v() const { return x_.template segment<D>(0); }
    VecD b() const { return x_.template segment<D>(D); }
    const Mat2D& P() const { return P_; }

private:
    Vec2D x_{Vec2D::Zero()};
    Mat2D P_{Mat2D::Identity()};
    MatD  R_{MatD::Identity() * 1e-3};
    double dt_k = 0.001;
    double scale_qv_{1.0};
    double scale_qb_{1e-4};
    double sigma_a2_{1e-4};
    static constexpr double QV_FLOOR = 1e-12;
    static constexpr double QB_FLOOR = 1e-12;
    double xddot_thresh = 4.0;
    double yddot_thresh = 6.0;
    double zddot_thresh = 7.0;
};
