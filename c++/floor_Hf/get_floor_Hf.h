#ifndef GET_FLOOR_HFH
#define GET_FLOOR_HFH

#include <Eigen/Dense>
#include "solver_floor_Hf.h"
#include "posedata.h"

PoseData get_floor_Hf(Eigen::MatrixXd &x1, Eigen::MatrixXd &x2, Eigen::Matrix3d &R1, Eigen::Matrix3d &q2, double f1);

#ifdef MATLAB_MEX_FILE /* This macro is defined automatically when using MATLAB */
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
#endif

double get_algebraic_error_floor_Hf(Eigen::VectorXd &data);

#endif /* GET_FLOOR_HFH */
