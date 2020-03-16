#include "get_floor_fHf.h"
#include "ransac_floor_fHf.h"
#include "ransac_data.h"
#include "randperm.h"
#include <Eigen/Geometry>
#include <vector>

using namespace Eigen;

RansacData ransac_floor_fHf(MatrixXd &p1, MatrixXd &p2, Matrix3d &R1, Matrix3d &R2, int nbr_iter, double thresh)
{
    // Init RANSAC loop
    PoseData best_posedata;
    int best_nbr_inliers = 0;
    VectorXi best_inliers;
    best_inliers.setZero();
    int nbr_pts = p1.cols();

    MatrixXd reproj1(2, nbr_pts);
    MatrixXd reproj2(2, nbr_pts);
    MatrixXd reproj(4, nbr_pts);
    VectorXd reproj_mean(nbr_pts);
    MatrixXd H(3,3);
    VectorXi inliers(nbr_pts);
    int nbr_inliers;
    VectorXi history(nbr_iter);

    int nbr_pts_minimal = 3;

    for (int i = 0; i < nbr_iter; i++) {
        std::vector<int> rands = randperm(nbr_pts_minimal, nbr_pts);
        MatrixXd x1(3, nbr_pts_minimal);
        MatrixXd x2(3, nbr_pts_minimal);
        for (int j = 0; j < nbr_pts_minimal; j++) {
            x1.col(j) << p1.col(rands[j]);
            x2.col(j) << p2.col(rands[j]);
        }

        // Get pose
        MatrixXd x1h(2,3);
        MatrixXd x2h(2,3);
        x1h << x1.colwise().hnormalized();
        x2h << x2.colwise().hnormalized();

        PoseData posedata = get_floor_fHf(x1h, x2h, R1, R2);

        // Compute reprojection error, and compare to other solutions.
        H = posedata.homography;
        reproj1 = p2.colwise().hnormalized() - (H * p1).colwise().hnormalized();
        reproj2 = p1.colwise().hnormalized() - (H.lu().solve(p2)).colwise().hnormalized();
        reproj << reproj1.cwiseAbs(), reproj2.cwiseAbs();
        reproj_mean = reproj.colwise().mean().array();

        inliers = (reproj_mean.array() < thresh).cast<int>();
        nbr_inliers = inliers.sum();

        if (best_nbr_inliers < nbr_inliers) {
            best_posedata = posedata;
            best_nbr_inliers = nbr_inliers;
            best_inliers = inliers;
        }
        history(i) = best_nbr_inliers;

    }

    RansacData ransac_data;
    ransac_data.posedata = best_posedata;
    ransac_data.inliers = best_inliers;
    ransac_data.history = history;

    return ransac_data;
}

// ---------------- //
// MATLAB interface //
// ---------------- //
#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("ransac_floor_fHf:nrhs", "Six input arguments are required.");
    }
    if (nlhs != 4) {
        mexErrMsgIdAndTxt("ransac_floor_fHf:nlhs", "Four output arguments are required.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("ransac_floor_fHf:notDouble", "Input data must be type double.");
    }
    if(mxGetNumberOfElements(prhs[0]) % 3 != 0
       && mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1])
       && mxGetNumberOfElements(prhs[2]) != 9 && mxGetNumberOfElements(prhs[3]) != 9
       && mxGetNumberOfElements(prhs[4]) != 1 && mxGetNumberOfElements(prhs[5]) != 1)
    {
        mexErrMsgIdAndTxt("ransac_floor_fHf:incorrectSize", "Input dimensions incorrect.");
    }
    // Convert to expected input
    int nbr_pts = mxGetNumberOfElements(prhs[0]) / 3;
    VectorXd x1_tmp = Map<VectorXd>(mxGetPr(prhs[0]), mxGetNumberOfElements(prhs[0]));
    VectorXd x2_tmp = Map<VectorXd>(mxGetPr(prhs[1]), mxGetNumberOfElements(prhs[1]));
    MatrixXd x1 = Map<MatrixXd>(x1_tmp.data(), 3, nbr_pts);
    MatrixXd x2 = Map<MatrixXd>(x2_tmp.data(), 3, nbr_pts);

    VectorXd R1_tmp = Map<VectorXd>(mxGetPr(prhs[2]), 9);
    VectorXd R2_tmp = Map<VectorXd>(mxGetPr(prhs[3]), 9);
    Matrix3d R1 = Map<Matrix3d>(R1_tmp.data(), 3, 3);
    Matrix3d R2 = Map<Matrix3d>(R2_tmp.data(), 3, 3);

    double *nbr_iter_p = mxGetPr(prhs[4]);
    double *thresh = mxGetPr(prhs[5]);
    int nbr_iter = (int) nbr_iter_p[0];

    // Compute output
    RansacData ransac_data = ransac_floor_fHf(x1, x2, R1, R2, nbr_iter, thresh[0]);
    PoseData posedata = ransac_data.posedata;

    // Wrap it up to Matlab compatible output
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* zr = mxGetPr(plhs[0]);
    for (Index i = 0; i < posedata.homography.size(); i++) {
        zr[i] = posedata.homography(i);
    }
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    zr = mxGetPr(plhs[1]);
    zr[0] = posedata.focal_length;

    plhs[2] = mxCreateDoubleMatrix(nbr_pts, 1, mxREAL);
    zr = mxGetPr(plhs[2]);
    for (Index i = 0; i < nbr_pts; i++) {
        zr[i] = ransac_data.inliers(i);
    }

    plhs[3] = mxCreateDoubleMatrix(nbr_iter, 1, mxREAL);
    zr = mxGetPr(plhs[3]);
    for (Index i = 0; i < nbr_iter; i++) {
        zr[i] = ransac_data.history(i);
    }
}
#endif
