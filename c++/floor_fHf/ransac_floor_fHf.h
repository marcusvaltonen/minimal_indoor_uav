#ifndef RANSAC_FLOOR_FHFH
#define RANSAC_FLOOR_FHFH

#include <Eigen/Dense>
#include "ransac_data.h"

RansacData ransac_floor_fHf(Eigen::MatrixXd &x1, Eigen::MatrixXd &x2, Eigen::Matrix3d &R1, Eigen::Matrix3d &q2, int nbr_iter, double threshold);

#ifdef MATLAB_MEX_FILE /* This macro is defined automatically when using MATLAB */
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
#endif

#endif /* RANSAC_FLOOR_FHFH */
