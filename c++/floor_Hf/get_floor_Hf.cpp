#include <float.h>
#include <cmath>
#include "get_floor_Hf.h"
#include "normalize2dpts.h"
#include <Eigen/Geometry>

using namespace Eigen;

PoseData get_floor_Hf(MatrixXd &p1, MatrixXd &p2, Matrix3d &R1, Matrix3d &R2, double f1)
{
    int nbr_coeffs = 21;
    int nbr_unknowns = 6;

    // Save copies of the inverse rotation
    Matrix3d R1T = R1.transpose();
    Matrix3d R2T = R2.transpose();

    // Compute normalization matrix
    double scale = normalize2dpts(p2);
    Vector3d s;
    s << scale, scale, 1.0;
    DiagonalMatrix<double, 3> S = s.asDiagonal();

    // Initialize known calibration matrix
    Vector3d k1inv;
    k1inv << 1 / f1, 1 / f1, 1.0;
    DiagonalMatrix<double, 3> K1inv = k1inv.asDiagonal();

    // Normalize data
    Matrix3d y1;
    Matrix3d x2;
    y1 = R1T * K1inv * p1.colwise().homogeneous();
    x2 = p2.colwise().homogeneous();

    x2 = S * x2;

    MatrixXd y1t(2,3);
    y1t << y1.colwise().hnormalized();
    MatrixXd x2t(2,3);
    x2t << x2.colwise().hnormalized();

    // Wrap input data to expected format
    VectorXd input(nbr_coeffs);
    input << y1t.col(0),
             x2t.col(0),
             y1t.col(1),
             x2t.col(1),
             y1t.col(2),
             x2t.col(2),
             Map<VectorXd>(R2T.data(), 9);

    // TODO: Not sure if this is necessary (assure const)
    const Map<VectorXd> input_data(input.data(), nbr_coeffs);

    // Extract solution
    MatrixXcd sols = solver_floor_Hf(input_data);

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
    ArrayXd xx(nbr_unknowns);
    VectorXd input_algebraic(nbr_coeffs + nbr_unknowns);

    for (int i = 0; i < real_sols.size(); i++) {
        if (real_sols(i) <= thresh) {
            // Compute algebraic error, and compare to other solutions.
            xx = sols.col(i).real();
            input_algebraic << xx, input;
            algebraic_error = get_algebraic_error_floor_Hf(input_algebraic);

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
    Matrix3d K, H;
    K = Vector3d(best_focal_length, best_focal_length, 1).asDiagonal();
    // Ki = Vector3d(1, 1, best_focal_length).asDiagonal();
    H = S.inverse() * K * R2 * best_homography * R1.transpose() * K1inv;

    // Package output
    PoseData posedata;
    posedata.homography = H;
    posedata.focal_length = best_focal_length / scale;

    return posedata;
}

// ---------------- //
// MATLAB interface //
// ---------------- //
#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 5) {
		mexErrMsgIdAndTxt("get_floor_Hf:nrhs", "Five input arguments are required.");
	}
	if (nlhs != 2) {
		mexErrMsgIdAndTxt("get_floor_Hf:nlhs", "Two output arguments are required.");
	}
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("get_floor_Hf:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) != 6 && mxGetNumberOfElements(prhs[1]) != 6 && mxGetNumberOfElements(prhs[2]) != 9 && mxGetNumberOfElements(prhs[3]) != 9 && mxGetNumberOfElements(prhs[4]) != 1) {
		mexErrMsgIdAndTxt("get_floor_Hf:incorrectSize", "Input dimensions incorrect.");
	}
    // Convert to expected input
    VectorXd x1_tmp = Map<VectorXd>(mxGetPr(prhs[0]), 6);
    VectorXd x2_tmp = Map<VectorXd>(mxGetPr(prhs[1]), 6);
    MatrixXd x1 = Map<MatrixXd>(x1_tmp.data(), 2, 3);
    MatrixXd x2 = Map<MatrixXd>(x2_tmp.data(), 2, 3);

    VectorXd R1_tmp = Map<VectorXd>(mxGetPr(prhs[2]), 9);
    VectorXd R2_tmp = Map<VectorXd>(mxGetPr(prhs[3]), 9);
    Matrix3d R1 = Map<Matrix3d>(R1_tmp.data(), 3, 3);
    Matrix3d R2 = Map<Matrix3d>(R2_tmp.data(), 3, 3);

    double *f1 = mxGetPr(prhs[4]);

    // Compute output
	PoseData posedata = get_floor_Hf(x1, x2, R1, R2, f1[0]);

    // Wrap it up to Matlab compatible output
	plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
	double* zr = mxGetPr(plhs[0]);
    for (Index i = 0; i < posedata.homography.size(); i++) {
        zr[i] = posedata.homography(i);
    }

    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    zr = mxGetPr(plhs[1]);
    zr[0] = posedata.focal_length;
}
#endif

// Function that utilizes the last equation of the DLT system to discard false solutions
double get_algebraic_error_floor_Hf(VectorXd &data)
{
    const double* d = data.data();

    // Compute algebraic error
    double error;
    error = -d[0]*d[5]*d[14]*d[25] - d[0]*d[14]*d[16]*d[19] - d[0]*d[14]*d[17]*d[22] - d[1]*d[5]*d[25] - d[1]*d[16]*d[19] - d[1]*d[17]*d[22] - d[2]*d[5]*d[15]*d[25] - d[2]*d[15]*d[16]*d[19] - d[2]*d[15]*d[17]*d[22] + d[3]*d[5]*d[15]*d[24] + d[3]*d[15]*d[16]*d[18] + d[3]*d[15]*d[17]*d[21];

    return abs(error);
}
