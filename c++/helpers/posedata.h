#ifndef POSEDATAH
#define POSEDATAH
#include <Eigen/Dense>
struct PoseData {
    Eigen::Matrix3d homography;
    double focal_length;
};

struct PoseDataVarFocal {
    Eigen::Matrix3d homography;
    double focal_length1;
    double focal_length2;
};
#endif /* POSEDATAH */
