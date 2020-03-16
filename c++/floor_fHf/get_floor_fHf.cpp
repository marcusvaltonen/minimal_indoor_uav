#include <float.h>
#include <cmath>
#include "get_floor_fHf.h"
#include "normalize2dpts.h"
#include <Eigen/Geometry>

using namespace Eigen;

PoseData get_floor_fHf(MatrixXd &p1, MatrixXd &p2, Matrix3d &R1, Matrix3d &R2)
{
    int nbr_coeffs = 30;
    int nbr_unknowns = 6;

    // Save copies of the inverse rotation
    Matrix3d R1T = R1.transpose();
    Matrix3d R2T = R2.transpose();

    // Compute normalization matrix
    double scale = normalize2dpts(p1);
    Vector3d s;
    s << scale, scale, 1.0;
    DiagonalMatrix<double, 3> S = s.asDiagonal();

    // Normalize data
    Matrix3d x1;
    Matrix3d x2;
    x1 = p1.colwise().homogeneous();
    x2 = p2.colwise().homogeneous();

    x1 = S * x1;
    x2 = S * x2;

    MatrixXd x1t(2,3);
    x1t << x1.colwise().hnormalized();
    MatrixXd x2t(2,3);
    x2t << x2.colwise().hnormalized();

    // Wrap input data to expected format
    VectorXd input(nbr_coeffs);
    input << x1t.col(0),
             x2t.col(0),
             x1t.col(1),
             x2t.col(1),
             x1t.col(2),
             x2t.col(2),
             Map<VectorXd>(R1T.data(), 9),
             Map<VectorXd>(R2T.data(), 9);

    // TODO: Not sure if this is necessary (assure const)
    const Map<VectorXd> input_data(input.data(), nbr_coeffs);

    // Extract solution
    MatrixXcd sols = solver_floor_fHf(input_data);

    // Pre-processing: Remove complex-valued solutions
    double thresh = 1e-5;
    ArrayXd real_sols(7);
    real_sols = sols.imag().cwiseAbs().colwise().sum();
    int nbr_real_sols = (real_sols <= thresh).count();

    // Allocate space for putative (real) homographies
    MatrixXd best_homography(3, 3);
    double best_focal_length;
    double best_algebraic_error = DBL_MAX;
    double algebraic_error;

    // Since this is a 2.5 pt solver, use the last
    // (previously unused) constraint, to discard
    // false solutions.
    ArrayXd xx(6);
    VectorXd input_algebraic(nbr_coeffs + nbr_unknowns);

    for (int i = 0; i < real_sols.size(); i++) {
        if (real_sols(i) <= thresh) {
            // Compute algebraic error, and compare to other solutions.
            xx = sols.col(i).real();
            input_algebraic << xx, input;
            algebraic_error = get_algebraic_error_floor_fHf(input_algebraic);

            if (algebraic_error < best_algebraic_error) {
                best_algebraic_error = algebraic_error;
                best_homography << xx[0], xx[2], xx[1],
                                       0, xx[3],     0,
                                  -xx[1], xx[4], xx[0];
                best_focal_length = xx[5];
            }
        }
    }
    // Construct homography
    Matrix3d K, Ki, H;
    K = Vector3d(best_focal_length, best_focal_length, 1).asDiagonal();
    Ki = Vector3d(1, 1, best_focal_length).asDiagonal();
    H = S.inverse() * K * R2 * best_homography * R1.transpose() * Ki * S;

    // Package output
    PoseData posedata;
    posedata.homography = H;
    posedata.focal_length = best_focal_length / scale;

    return posedata;
}

// Function that utilizes the last equation of the DLT system to discard false solutions
double get_algebraic_error_floor_fHf(VectorXd &data)
{
    const double* d = data.data();

    // Compute algebraic error
    double error;
    error = -d[0]*std::pow(d[5],2)*d[24]*d[34] - d[0]*d[5]*d[14]*d[18]*d[34] - d[0]*d[5]*d[15]*d[21]*d[34] - d[0]*d[5]*d[16]*d[24]*d[28] - d[0]*d[5]*d[17]*d[24]*d[31] - d[0]*d[14]*d[16]*d[18]*d[28] - d[0]*d[14]*d[17]*d[18]*d[31] - d[0]*d[15]*d[16]*d[21]*d[28] - d[0]*d[15]*d[17]*d[21]*d[31] - d[1]*std::pow(d[5],2)*d[26]*d[34] - d[1]*d[5]*d[14]*d[20]*d[34] - d[1]*d[5]*d[15]*d[23]*d[34] - d[1]*d[5]*d[16]*d[26]*d[28] - d[1]*d[5]*d[17]*d[26]*d[31] - d[1]*d[14]*d[16]*d[20]*d[28] - d[1]*d[14]*d[17]*d[20]*d[31] - d[1]*d[15]*d[16]*d[23]*d[28] - d[1]*d[15]*d[17]*d[23]*d[31] - d[2]*std::pow(d[5],2)*d[25]*d[34] - d[2]*d[5]*d[14]*d[19]*d[34] - d[2]*d[5]*d[15]*d[22]*d[34] - d[2]*d[5]*d[16]*d[25]*d[28] - d[2]*d[5]*d[17]*d[25]*d[31] - d[2]*d[14]*d[16]*d[19]*d[28] - d[2]*d[14]*d[17]*d[19]*d[31] - d[2]*d[15]*d[16]*d[22]*d[28] - d[2]*d[15]*d[17]*d[22]*d[31] + d[3]*std::pow(d[5],2)*d[25]*d[33] + d[3]*d[5]*d[14]*d[19]*d[33] + d[3]*d[5]*d[15]*d[22]*d[33] + d[3]*d[5]*d[16]*d[25]*d[27] + d[3]*d[5]*d[17]*d[25]*d[30] + d[3]*d[14]*d[16]*d[19]*d[27] + d[3]*d[14]*d[17]*d[19]*d[30] + d[3]*d[15]*d[16]*d[22]*d[27] + d[3]*d[15]*d[17]*d[22]*d[30];
    return abs(error);
}
