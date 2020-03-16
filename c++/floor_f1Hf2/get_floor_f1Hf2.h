#ifndef GET_FLOOR_F1HF2H
#define GET_FLOOR_F1HF2H

#include <vector>
#include <Eigen/Dense>
#include "solver_floor_f1Hf2.h"
#include "posedata.h"

std::vector<PoseDataVarFocal> get_floor_f1Hf2(Eigen::MatrixXd &x1, Eigen::MatrixXd &x2, Eigen::Matrix3d &R1, Eigen::Matrix3d &R2);

#ifdef MATLAB_MEX_FILE /* This macro is defined automatically when using MATLAB */
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
#endif

#endif /* GET_FLOOR_F1HF2H */
