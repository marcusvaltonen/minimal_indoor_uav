#include <iostream>
#include <float.h>
#include <cmath>
#include "get_floor_f1Hf2.h"
#include "normalize2dpts.h"
#include <Eigen/Geometry>

using namespace Eigen;

std::vector<PoseDataVarFocal> get_floor_f1Hf2(MatrixXd &p1, MatrixXd &p2, Matrix3d &R1, Matrix3d &R2)
{
    int nbr_coeffs = 9 * 5;
    int nbr_unknowns = 4;

    // Save copies of the inverse rotation
    Matrix3d R1T = R1.transpose();
    Matrix3d R2T = R2.transpose();

    // Compute normalization matrices
    double scale1 = normalize2dpts(p1);
    double scale2 = normalize2dpts(p2);
    Vector3d s1, s2;
    s1 << scale1, scale1, 1.0;
    s2 << scale2, scale2, 1.0;
    DiagonalMatrix<double, 3> S1 = s1.asDiagonal();
    DiagonalMatrix<double, 3> S2 = s2.asDiagonal();

    // Normalize data
    MatrixXd x1(3, 3);
    MatrixXd x2(3, 3);
    x1 = S1 * p1;
    x2 = S2 * p2;

    // Setup DLT equations
    MatrixXd A(6,9);
    A << 0, 0, 0, -x2(2,0)*x1.col(0).transpose(), x2(1,0)*x1.col(0).transpose(),
         x2(2,0)*x1.col(0).transpose(), 0, 0, 0, -x2(0,0)*x1.col(0).transpose(),
         0, 0, 0, -x2(2,1)*x1.col(1).transpose(), x2(1,1)*x1.col(1).transpose(),
         x2(2,1)*x1.col(1).transpose(), 0, 0, 0, -x2(0,1)*x1.col(1).transpose(),
         0, 0, 0, -x2(2,2)*x1.col(2).transpose(), x2(1,2)*x1.col(2).transpose(),
         x2(2,2)*x1.col(2).transpose(), 0, 0, 0, -x2(0,2)*x1.col(2).transpose();

    JacobiSVD<MatrixXd> svd(A, ComputeFullV);
    ArrayXXd V = svd.matrixV();

    // Wrap input data to expected format
    VectorXd input(nbr_coeffs);
    input << V.col(6),
             V.col(7),
             V.col(8),
             Map<VectorXd>(R1T.data(), 9),
             Map<VectorXd>(R2T.data(), 9);

    // TODO: Not sure if this is necessary (assure const)
    const Map<VectorXd> input_data(input.data(), nbr_coeffs);

    // Extract solution
    MatrixXcd sols = solver_floor_f1Hf2(input_data);

    // Pre-processing: Remove complex-valued solutions
    double thresh = 1e-5;
    ArrayXd real_sols(5);
    real_sols = sols.imag().cwiseAbs().colwise().sum();
    int nbr_real_sols = (real_sols <= thresh).count();

    // Allocate space for putative (real) homographies

    // Since this is a 4 pt solver, we only return the solutions.
    std::vector<PoseDataVarFocal> posedata(nbr_real_sols);
    ArrayXd xx(nbr_unknowns);
    Matrix3d Htmp;
    double f1, f2;
    int cnt = 0;

    for (int i = 0; i < real_sols.size(); i++) {
        if (real_sols(i) <= thresh) {
            xx = sols.col(i).real();

            // Extract focal lengths
            f1 = 1 / xx[2] / scale1;
            f2 = xx[3] / scale2;

            // Construct putative homography
            VectorXd tmp(9);
            tmp << V.col(6) + xx[0] * V.col(7) + xx[1] * V.col(8);
            Htmp = Map<Matrix3d>(tmp.data(), 3 ,3);
            Htmp.transposeInPlace();
            Htmp = S2.inverse() * Htmp * S1;

            // Append
            posedata[cnt].focal_length1 = f1;
            posedata[cnt].focal_length2 = f2;
            posedata[cnt].homography = Htmp;
            cnt++;
        }
    }

    return posedata;
}

// ---------------- //
// MATLAB interface //
// ---------------- //

#ifdef MATLAB_MEX_FILE /* This macro is defined automatically when using MATLAB */
#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 4) {
		mexErrMsgIdAndTxt("get_floor_f1Hf2:nrhs", "Four input arguments are required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("get_floor_f1Hf2:nlhs", "One output arguments is required.");
	}
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("get_floor_f1Hf2:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) != 9 && mxGetNumberOfElements(prhs[1]) != 9 && mxGetNumberOfElements(prhs[2]) != 9 && mxGetNumberOfElements(prhs[3]) != 9) {
		mexErrMsgIdAndTxt("get_floor_f1Hf2:incorrectSize", "Input dimensions incorrect.");
	}
    // Convert to expected input
    VectorXd x1_tmp = Map<VectorXd>(mxGetPr(prhs[0]), 9);
    VectorXd x2_tmp = Map<VectorXd>(mxGetPr(prhs[1]), 9);
    MatrixXd x1 = Map<MatrixXd>(x1_tmp.data(), 3, 3);
    MatrixXd x2 = Map<MatrixXd>(x2_tmp.data(), 3, 3);

    VectorXd R1_tmp = Map<VectorXd>(mxGetPr(prhs[2]), 9);
    VectorXd R2_tmp = Map<VectorXd>(mxGetPr(prhs[3]), 9);
    Matrix3d R1 = Map<Matrix3d>(R1_tmp.data(), 3, 3);
    Matrix3d R2 = Map<Matrix3d>(R2_tmp.data(), 3, 3);

    // Compute output
	std::vector<PoseDataVarFocal> posedata = get_floor_f1Hf2(x1, x2, R1, R2);

    // Wrap it all up
    std::size_t NUMBER_OF_STRUCTS = posedata.size();
    const char *field_names[] = {"H", "f1", "f2"};
    mwSize dims[2] = {1, NUMBER_OF_STRUCTS };
    int H_field, f1_field, f2_field;
    mwIndex i;

    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);

    H_field = mxGetFieldNumber(plhs[0], "H");
    f1_field = mxGetFieldNumber(plhs[0], "f1");
    f2_field = mxGetFieldNumber(plhs[0], "f2");

	double* zr;
    for (i = 0; i < NUMBER_OF_STRUCTS; i++) {
        mxArray *field_value;

        // Create H
        field_value = mxCreateDoubleMatrix(3, 3, mxREAL);
        zr = mxGetPr(field_value);
        for (Index j = 0; j < posedata[i].homography.size(); j++) {
            zr[j] = posedata[i].homography(j);
        }
        mxSetFieldByNumber(plhs[0],i,H_field,field_value);

        // Create f1
        field_value = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
        zr = mxGetPr(field_value);
        zr[0] = posedata[i].focal_length1;
        mxSetFieldByNumber(plhs[0],i,f1_field,field_value);

        // Create f2
        field_value = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
        zr = mxGetPr(field_value);
        zr[0] = posedata[i].focal_length2;
        mxSetFieldByNumber(plhs[0],i,f2_field,field_value);
    }
}
#endif
