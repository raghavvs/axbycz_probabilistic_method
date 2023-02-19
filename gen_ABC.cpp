#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <random>

using namespace Eigen;
using namespace std;

// Generate a random transformation matrix from a 6D vector
Matrix4d expm(const VectorXd& t) {
    Matrix4d T;
    T.setIdentity();
    double theta = t.head<3>().norm();
    if (theta < 1e-12) {
        // Use first-order Taylor expansion for small angles
        T.block<3, 1>(0, 3) = t.tail<3>();
        T.block<3, 3>(0, 0) = Matrix3d::Identity() + skew(t.tail<3>());
    } else {
        Vector3d axis = t.head<3>() / theta;
        T.block<3, 3>(0, 0) = AngleAxisd(theta, axis).toRotationMatrix();
        T.block<3, 1>(0, 3) = (Matrix3d::Identity() - T.block<3, 3>(0, 0)) * skew(axis) * t.tail<3>();
    }
    return T;
}

// Skew-symmetric matrix from a 3D vector
Matrix3d skew(const Vector3d& v) {
    Matrix3d S;
    S << 0, -v(2), v(1),
         v(2), 0, -v(0),
         -v(1), v(0), 0;
    return S;
}

// Generate a random 6D vector from a multivariate Gaussian distribution
VectorXd mvg(const MatrixXd& M, const MatrixXd& Sig) {
    std::random_device rd;
    std::mt19937 gen(rd());
    MatrixXd R = M + Sig.llt().matrixL() * MatrixXd::NullaryExpr(Sig.cols(), 1, [&]() { return std::normal_distribution<>(0.0)(gen); });
    return R;
}

void ABC_Generate(int len, int opt, const Matrix4d& A_initial, const Matrix4d& B_initial, const Matrix4d& C_initial, const Matrix4d& X, const Matrix4d& Y, const Matrix4d& Z, Matrix4d& A, Matrix4d& B, Matrix4d& C) {
    Matrix4d I4 = Matrix4d::Identity();
    Matrix4d Y_inv = Y.inverse();
    Matrix4d X_inv = X.inverse();
    Matrix4d Z_inv = Z.inverse();
    MatrixXd M(6, 1);
    M.setZero();
    MatrixXd Sig(6, 6);
    Sig.setIdentity();

    if (opt == 1) { // Fix A, randomize B and C
        B.resize(4, 4 * len);
        C.resize(4, 4 * len);
        for (int m = 0; m < len; m++) {
            Matrix4d B_m = expm(mvg(M, Sig)) * B_initial;
            Matrix4d C_m = Y_inv * A_initial * X * B_m * Z_inv;
            B.block<4, 4>(0, 4 * m) = B_m;
            C.block<4, 4>(0, 4 * m) = C_m;
        }
        A = A_initial;
    } else if (opt == 2) { // Fix B, randomize A and C
        A.resize(4, 4 * len);
        C.resize(4,


// Helper function to convert 6x1 se(3) vector to 4x4 homogeneous transformation matrix
Matrix4d expm_se3(Vector6d v) {
    double theta = v.head(3).norm();
    Matrix3d w_hat;
    w_hat << 0, -v(2), v(1),
             v(2), 0, -v(0),
            -v(1), v(0), 0;
    Matrix3d R;
    if (theta < 1e-10) {
        R = Matrix3d::Identity() + w_hat;
    } else {
        R = Matrix3d::Identity() + sin(theta)/theta * w_hat + (1-cos(theta))/(theta*theta) * w_hat*w_hat;
    }
    Vector3d p = (Matrix3d::Identity() - R) * v.tail(3);
    Matrix4d T = Matrix4d::Identity();
    T.block(0, 0, 3, 3) = R;
    T.block(0, 3, 3, 1) = p;
    return T;
}

// Helper function to generate random 6x1 vector from a multivariate Gaussian distribution
Vector6d mvg(Vector3d mu, Matrix3d Sigma) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::multivariate_normal_distribution<double> dist(mu.data(), Sigma.data());
    Vector6d v;
    for (int i = 0; i < 6; i++) {
        v(i) = dist(gen);
    }
    return v;
}

// Data generation for AXB = YCZ problem
void ABC_Generate(int length, int opt, Vector3d M, Matrix3d Sig, Matrix4d X, Matrix4d Y, Matrix4d Z, Matrix4d& A, Matrix4d& B, Matrix4d& C) {
    int len = length; // number of generated data pairs
    Matrix4d A_initial, B_initial, C_initial;

    // Generate initial transformations A, B, and C
    std::vector<Vector3d> qz1 { pi/6, pi/3, pi/4, pi/4, -pi/4, 0 };
    std::vector<Vector3d> qz2 { pi/3, pi/4, pi/3, -pi/4, pi/4, 0 };
    std::vector<Vector3d> qz3 { pi/4, pi/3, -pi/3, -pi/6, pi/4, 0 };

    Vector6d a, b, c;
    std::normal_distribution<double> norm_dist(0.0, 1.0);

    for (int i = 0; i < 6; i++) {
        a(i) = norm_dist(gen);
        b(i) = norm_dist(gen);
        c(i) = norm_dist(gen);
    }

    a /= a.norm();
    b /= b.norm();
    c /= c.norm();

    A_initial = expm_se3(a);
    B_initial = expm_se3(b);
    C_initial = expm_se3(c);

    if (opt == 1) {
        // Fix A, randomize B and C
        B.resize(4, 4, len);
        C.resize(4, 4, len);

        for (int m = 0


Matrix4d hat(Matrix<double, 6, 1> omega)
{
    Matrix4d omega_hat;
    omega_hat << 0, -omega(2), omega(1), omega(3),
                 omega(2), 0, -omega(0), omega(4),
                 -omega(1), omega(0), 0, omega(5),
                 0, 0, 0, 0;
    return omega_hat;
}

Matrix4d exp_hat(Matrix<double, 6, 1> omega, double theta)
{
    Matrix<double, 6, 1> omega_norm = omega / omega.norm();
    Matrix4d omega_hat = hat(omega_norm);
    Matrix4d R = Matrix4d::Identity() + omega_hat * sin(theta) + omega_hat * omega_hat * (1 - cos(theta));
    return R;
}

Matrix<double, 6, 1> log_hat(Matrix4d R)
{
    double trace_R = R(0,0) + R(1,1) + R(2,2);
    double theta = acos((trace_R - 1) / 2);

    Matrix<double, 6, 1> omega;
    if (abs(theta) < 1e-10) {
        omega << 0, 0, 0, R(0,3), R(1,3), R(2,3);
    } else if (abs(M_PI - theta) < 1e-10) {
        Matrix4d A = R - Matrix4d::Identity();
        Matrix<double, 3, 1> v;
        if (abs(A(0,0)) > 1e-10 && abs(A(1,1)) > 1e-10) {
            v << A(0,2), A(1,2), A(1,0);
        } else if (abs(A(0,0)) > 1e-10 && abs(A(2,2)) > 1e-10) {
            v << A(2,1), A(0,2), A(2,0);
        } else {
            v << A(1,0), A(2,1), A(0,1);
        }
        double s = sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2));
        if (s > 1e-10) {
            v /= s;
        }
        omega << v(0), v(1), v(2), R(0,3), R(1,3), R(2,3);
    } else {
        Matrix4d A = (R - R.transpose()) / (2 * sin(theta));
        omega << A(2,1), A(0,2), A(1,0), R(0,3), R(1,3), R(2,3);
    }
    return theta * omega / omega.norm();
}

Matrix4d ABC_Generate(int len, int opt, Matrix<double, 6, 1> M, Matrix<double, 6, 6> Sig, Matrix4d X, Matrix4d Y, Matrix4d Z)
{
    Matrix<double, 6, 1> a, b, c;
    Matrix4d A_initial, B_initial, C_initial;

    //using puma560 to generate tranformation A


// Function to generate random 4x4 transformation matrix using axis-angle representation
Matrix4d expm_se3(const VectorXd& omega)
{
    Matrix3d Omega;
    Omega << 0, -omega(2), omega(1), omega(2), 0, -omega(0), -omega(1), omega(0), 0;
    double theta = omega.norm();
    Matrix3d Omega_hat;
    if (theta < 1e-10) {
        Omega_hat.setZero();
    }
    else {
        Omega_hat = Omega / theta;
    }
    Matrix4d T;
    T.topLeftCorner(3, 3) = Matrix3d::Identity() + std::sin(theta) * Omega_hat + (1 - std::cos(theta)) * Omega_hat * Omega_hat;
    T.topRightCorner(3, 1) = (Matrix3d::Identity() - T.topLeftCorner(3, 3)) * omega.cross(omega.tail<3>()) + Omega_hat * omega.tail<3>() * omega.transpose();
    T.bottomLeftCorner(1, 3).setZero();
    T(3, 3) = 1;
    return T;
}

// Function to generate random matrix with Gaussian noise
VectorXd mvg(const VectorXd& mu, const MatrixXd& sigma, std::default_random_engine& generator)
{
    int n = mu.size();
    std::normal_distribution<double> distribution(0.0, 1.0);
    VectorXd x(n);
    for (int i = 0; i < n; i++) {
        x(i) = distribution(generator);
    }
    return mu + sigma.llt().matrixL() * x;
}

// Function to generate data for AXB = YCZ problem
void ABC_Generate(int len, int opt, const Matrix4d& M, const MatrixXd& Sig, const Matrix4d& X, const Matrix4d& Y, const Matrix4d& Z, Matrix4d& A, Matrix4d& B, Matrix4d& C)
{
    std::default_random_engine generator;
    generator.seed(0);

    Matrix4d A_initial = Matrix4d::Identity();
    Matrix4d B_initial = Matrix4d::Identity();
    Matrix4d C_initial = Matrix4d::Identity();
    VectorXd a(6), b(6), c(6);
    a.setRandom();
    b.setRandom();
    c.setRandom();
    a /= a.norm();
    b /= b.norm();
    c /= c.norm();
    A_initial = expm_se3(a);
    B_initial = expm_se3(b);
    C_initial = expm_se3(c);

    if (opt == 1) { // Fix A, randomize B and C
        B = Matrix4d::Zero();
        C = Matrix4d::Zero();
        for (int m = 0; m < len; m++) {
            VectorXd b_noise = mvg(VectorXd::Zero(6), Sig, generator);
            Matrix4d B_temp = expm_se3(b_noise) * B_initial;
            C.col(m) = Y.inverse() * A_initial * X * B_temp * Z.inverse();
            B.col(m) = B_temp.col(3);
        }
        A = A_initial;
    }
    else if (opt == 2) { // Fix B, randomize A and C
        A = Matrix4d::Zero();
        C = Matrix4d::

// Generate a random transformation matrix from a 6D vector
Eigen::Matrix4d expm(se3_vector w) {
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    Eigen::Matrix3d W = skew(w.head(3));
    double theta = w.tail(3).norm();
    if (theta < 1e-12) {
        T.block<3, 3>(0, 0) += W;
        T.block<3, 1>(0, 3) = w.head(3);
    } else {
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + (sin(theta) / theta) * W + ((1 - cos(theta)) / (theta * theta)) * W * W;
        Eigen::Vector3d p = ((Eigen::Matrix3d::Identity() - R) * skew(w.head(3)) + w.tail(3) * w.head(3).transpose()) / (theta * theta);
        T.block<3, 3>(0, 0) = R;
        T.block<3, 1>(0, 3) = p;
    }
    return T;
}

// Generate a random transformation matrix using Puma560
Eigen::Matrix4d generate_puma560(int idx) {
    if (idx == 0) {
        double qz[] = {M_PI / 6, M_PI / 3, M_PI / 4, M_PI / 4, -M_PI / 4, 0};
        return p560.fkine(qz);
    } else if (idx == 1) {
        double qz[] = {M_PI / 3, M_PI / 4, M_PI / 3, -M_PI / 4, M_PI / 4, 0};
        return p560.fkine(qz);
    } else {
        double qz[] = {M_PI / 4, M_PI / 3, -M_PI / 3, -M_PI / 6, M_PI / 4, 0};
        return p560.fkine(qz);
    }
}

void ABC_Generate(int length, int opt, const Eigen::VectorXd& M, const Eigen::MatrixXd& Sig, const Eigen::Matrix4d& X, const Eigen::Matrix4d& Y, const Eigen::Matrix4d& Z, Eigen::Matrix4d& A, Eigen::Matrix4d& B, Eigen::Matrix4d& C) {
    // Generate a random transformation matrix A, B, and C
    Eigen::Vector3d a = Eigen::Vector3d::Random();
    a.normalize();
    A = expm(se3_vector(a));
    Eigen::Vector3d b = Eigen::Vector3d::Random();
    b.normalize();
    B = expm(se3_vector(b));
    Eigen::Vector3d c = Eigen::Vector3d::Random();
    c.normalize();
    C = expm(se3_vector(c));

    if (opt == 1) { // Fix A, randomize B and C
        B.resize(4, 4 * length);
        C.resize(4, 4 * length);

        for (int m = 0; m < length; m++) {
            Eigen::Vector3d bm = mvg(M, Sig);
            B.block<4, 4>(0, 4 * m) = expm(se3_vector(bm)) * B;
            C.block<4, 4>(0, 4 * m) = Y

        // calculate B using X, Y, Z, and A
        for (int m = 0; m < len; m++) {
            auto temp = mvg(M, Sig, 1);
            Mat se3_vec_temp = se3_vec(temp);

            // calculate B
            Mat expm_temp = expm(se3_vec_temp);
            Mat B_temp = expm_temp * B_initial;

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    Mat temp1 = Y.col(i);
                    Mat temp2 = A_initial * X.col(j);
                    Mat temp3 = B_temp.row(j);
                    Mat temp4 = Z.col(i);
                    B(i, j, m) = temp1.dot(temp2.dot(temp3.dot(temp4)));
                }
            }
        }

        A = A_initial;

    } else if (opt == 2) { // Fix B, randomize A and C
        // This can be applied to both serial-parallel and dual-robot arm calibrations

        Mat A(len, 4, 4);
        Mat C(len, 4, 4);

        // calculate A and C using B
        for (int m = 0; m < len; m++) {
            auto temp = mvg(M, Sig, 1);
            Mat se3_vec_temp = se3_vec(temp);

            // calculate A
            Mat expm_temp = expm(se3_vec_temp);
            A.row(m) = expm_temp * A_initial;

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    Mat temp1 = Y.col(i);
                    Mat temp2 = A(m).col(j);
                    Mat temp3 = X.col(j);
                    Mat temp4 = B_initial.row(j);
                    Mat temp5 = Z.col(i);
                    C(i, j, m) = temp1.dot(temp2.dot(temp3.dot(temp4.dot(temp5))));
                }
            }
        }

        B = B_initial;

    } else if (opt == 3) { // Fix C, randomize A and B
        // This is only physically achievable on multi-robot hand-eye calibration

        Mat A(len, 4, 4);
        Mat B(len, 4, 4);

        // calculate A and B using C
        for (int m = 0; m < len; m++) {
            auto temp = mvg(M, Sig, 1);
            Mat se3_vec_temp = se3_vec(temp);

            // calculate A
            Mat expm_temp = expm(se3_vec_temp);
            A.row(m) = expm_temp * A_initial;

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    Mat temp1 = X.col(j);
                    Mat temp2 = A(m).col(j);
                    Mat temp3 = Y.col(i);
                    Mat temp4 = C_initial.row(i);
                    Mat temp5 = Z.col(i);
                    B(i, j, m) = temp1.dot(temp2.dot(temp3.dot(temp4.dot(temp5))));
                }
            }
        }

        C = C_initial;

    } else if (opt == 4) { // This is for testing traditional AXBYCZ solver that demands the correspondence between the data pairs {A_i, B_i, C_i}

        Mat A(len,
        if (opt == 3) { // Fix C, randomize A and B
            A.resize(len);
            B.resize(len);
            C = {C_initial};
            for (int m = 0; m < len; m++) {
                A[m] = expm(se3_vec(mvg(M, Sig, 1))) * A_initial;
                // zg: need to have brackets when using mldivide notation
                B[m] = X.inv() * A[m].inv() * Y * C[0] * Z;
            }
        } else if (opt == 4) { // For testing traditional AXBYCZ solver that demands the correspondence between the data pairs {A_i, B_i, C_i}
            A.resize(len);
            B.resize(len);
            C.resize(len);
            for (int m = 0; m < len; m++) {
                A[m] = expm(se3_vec(mvg(M, Sig, 1))) * A_initial;
                C[m] = expm(se3_vec(mvg(M, Sig, 1))) * C_initial;
                // zg: need to have brackets when using mldivide notation
                B[m] = X.inv() * A[m].inv() * (Y * C[m] * Z);
                // B[m] = X.inv() * A[m].inv() * Y * C[m] * Z.inv();
            }
        }
        return std::make_tuple(A, B, C);
    }
}
