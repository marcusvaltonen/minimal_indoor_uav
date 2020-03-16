#ifndef GET_FLOOR_FHFH
#define GET_FLOOR_FHFH

#include <Eigen/Dense>
#include "solver_floor_fHf.h"
#include "posedata.h"

PoseData get_floor_fHf(Eigen::MatrixXd &x1, Eigen::MatrixXd &x2, Eigen::Matrix3d &R1, Eigen::Matrix3d &q2);

double get_algebraic_error_floor_fHf(Eigen::VectorXd &data);

#endif /* GET_FLOOR_FHFH */
