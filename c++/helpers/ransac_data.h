#ifndef RANSAC_DATAH
#define RANSAC_DATAH
#include <Eigen/Dense>
#include <vector>
#include <posedata.h>

struct RansacData {
    PoseData posedata;
    Eigen::VectorXi inliers;
    Eigen::VectorXi history;
};

#endif /* RANSAC_DATAH */
