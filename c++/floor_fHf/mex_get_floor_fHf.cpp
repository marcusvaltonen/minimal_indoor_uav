#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "posedata.h"
#include "get_floor_fHf.h"
#include <Eigen/Dense>

using namespace Eigen;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 4) {
		mexErrMsgIdAndTxt("get_floor_fHf:nrhs", "Four input arguments are required.");
	}
	if (nlhs != 2) {
		mexErrMsgIdAndTxt("get_floor_fHf:nlhs", "Two output arguments are required.");
	}
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("get_floor_fHf:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) != 6 && mxGetNumberOfElements(prhs[1]) != 6 && mxGetNumberOfElements(prhs[2]) != 9 && mxGetNumberOfElements(prhs[3]) != 9) {
		mexErrMsgIdAndTxt("get_floor_fHf:incorrectSize", "Input dimensions incorrect.");
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

    // Compute output
	PoseData posedata = get_floor_fHf(x1, x2, R1, R2);

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

