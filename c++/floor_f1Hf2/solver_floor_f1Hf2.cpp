#include <Eigen/Dense>

using namespace Eigen;

MatrixXcd solver_floor_f1Hf2(const VectorXd& data)
{
    // Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(48);

    double t2 = d[11]*d[37];
    double t3 = d[14]*d[40];
    double t4 = t2+t3;
    double t5 = d[20]*d[37];
    double t6 = d[23]*d[40];
    double t7 = t5+t6;
    double t8 = d[2]*d[37];
    double t9 = d[5]*d[40];
    double t10 = t8+t9;
    double t11 = d[33]*d[44];
    double t12 = d[35]*d[42];
    double t13 = t11+t12;
    double t14 = d[33]*d[42];
    double t16 = d[35]*d[44];
    double t15 = t14-t16;

    coeffs <<
    d[17]*d[35]*d[43],
    d[26]*d[35]*d[43],
    t4*d[35],
    t7*d[35],
    d[43]*(d[15]*d[29]+d[16]*d[32]),
    d[43]*(d[24]*d[29]+d[25]*d[32]),
    d[8]*d[35]*d[43],
    d[9]*d[29]*d[37]+d[10]*d[32]*d[37]+d[12]*d[29]*d[40]+d[13]*d[32]*d[40],
    d[18]*d[29]*d[37]+d[19]*d[32]*d[37]+d[21]*d[29]*d[40]+d[22]*d[32]*d[40],
    t10*d[35],
    d[43]*(d[6]*d[29]+d[7]*d[32]),
    d[0]*d[29]*d[37]+d[1]*d[32]*d[37]+d[3]*d[29]*d[40]+d[4]*d[32]*d[40],
    d[17]*d[33]*d[43],
    d[26]*d[33]*d[43],
    t4*d[33],
    t7*d[33],
    d[43]*(d[15]*d[27]+d[16]*d[30]),
    d[43]*(d[24]*d[27]+d[25]*d[30]),
    d[8]*d[33]*d[43],
    d[9]*d[27]*d[37]+d[10]*d[30]*d[37]+d[12]*d[27]*d[40]+d[13]*d[30]*d[40],
    d[18]*d[27]*d[37]+d[19]*d[30]*d[37]+d[21]*d[27]*d[40]+d[22]*d[30]*d[40],
    t10*d[33],
    d[43]*(d[6]*d[27]+d[7]*d[30]),
    d[0]*d[27]*d[37]+d[1]*d[30]*d[37]+d[3]*d[27]*d[40]+d[4]*d[30]*d[40],
    t13*d[17],
    t13*d[26],
    d[11]*d[33]*d[38]+d[11]*d[35]*d[36]+d[14]*d[33]*d[41]+d[14]*d[35]*d[39],
    d[20]*d[33]*d[38]+d[20]*d[35]*d[36]+d[23]*d[33]*d[41]+d[23]*d[35]*d[39],
    d[15]*d[27]*d[44]+d[15]*d[29]*d[42]+d[16]*d[30]*d[44]+d[16]*d[32]*d[42],
    d[24]*d[27]*d[44]+d[24]*d[29]*d[42]+d[25]*d[30]*d[44]+d[25]*d[32]*d[42],
    t13*d[8],
    d[9]*d[27]*d[38]+d[9]*d[29]*d[36]+d[10]*d[30]*d[38]+d[10]*d[32]*d[36]+d[12]*d[27]*d[41]+d[12]*d[29]*d[39]+d[13]*d[30]*d[41]+d[13]*d[32]*d[39],
    d[18]*d[27]*d[38]+d[18]*d[29]*d[36]+d[19]*d[30]*d[38]+d[19]*d[32]*d[36]+d[21]*d[27]*d[41]+d[21]*d[29]*d[39]+d[22]*d[30]*d[41]+d[22]*d[32]*d[39],
    d[2]*d[33]*d[38]+d[2]*d[35]*d[36]+d[5]*d[33]*d[41]+d[5]*d[35]*d[39],
    d[6]*d[27]*d[44]+d[6]*d[29]*d[42]+d[7]*d[30]*d[44]+d[7]*d[32]*d[42],
    d[0]*d[27]*d[38]+d[0]*d[29]*d[36]+d[1]*d[30]*d[38]+d[1]*d[32]*d[36]+d[3]*d[27]*d[41]+d[3]*d[29]*d[39]+d[4]*d[30]*d[41]+d[4]*d[32]*d[39],
    t15*d[17],
    t15*d[26],
    d[11]*d[33]*d[36]-d[11]*d[35]*d[38]+d[14]*d[33]*d[39]-d[14]*d[35]*d[41],
    d[20]*d[33]*d[36]-d[20]*d[35]*d[38]+d[23]*d[33]*d[39]-d[23]*d[35]*d[41],
    d[15]*d[27]*d[42]-d[15]*d[29]*d[44]+d[16]*d[30]*d[42]-d[16]*d[32]*d[44],
    d[24]*d[27]*d[42]-d[24]*d[29]*d[44]+d[25]*d[30]*d[42]-d[25]*d[32]*d[44],
    t15*d[8],
    d[9]*d[27]*d[36]-d[9]*d[29]*d[38]+d[10]*d[30]*d[36]+d[12]*d[27]*d[39]-d[10]*d[32]*d[38]-d[12]*d[29]*d[41]+d[13]*d[30]*d[39]-d[13]*d[32]*d[41],
    d[18]*d[27]*d[36]-d[18]*d[29]*d[38]+d[19]*d[30]*d[36]+d[21]*d[27]*d[39]-d[19]*d[32]*d[38]-d[21]*d[29]*d[41]+d[22]*d[30]*d[39]-d[22]*d[32]*d[41],
    d[2]*d[33]*d[36]-d[2]*d[35]*d[38]+d[5]*d[33]*d[39]-d[5]*d[35]*d[41],
    d[6]*d[27]*d[42]-d[6]*d[29]*d[44]+d[7]*d[30]*d[42]-d[7]*d[32]*d[44],
    d[0]*d[27]*d[36]-d[0]*d[29]*d[38]+d[1]*d[30]*d[36]+d[3]*d[27]*d[39]-d[1]*d[32]*d[38]-d[3]*d[29]*d[41]+d[4]*d[30]*d[39]-d[4]*d[32]*d[41];



    // Setup elimination template
    static const int coeffs0_ind[] = { 0,12,0,12,24,36,4,2,0,14,12,16,24,26,36,38,4,16,28,40,7,4,19,16,28,31,40,43,5,3,1,15,13,17,25,27,37,39,7,2,14,19,26,38,8,3,15,20,27,39,10,9,6,21,18,22,30,33,42,45,11,9,21,23,33,45 };
    static const int coeffs1_ind[] = { 8,5,20,17,29,32,41,44,7,19,31,43,8,20,32,44,11,10,23,22,34,35,46,47,11,23,35,47 };
    static const int C0_ind[] = { 0,5,11,13,17,19,20,21,22,23,24,25,26,27,28,29,31,33,37,39,41,42,43,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,62,64,65,66,68,70,72,74,75,76,78,80,81,82,83,84,85,86,87,88,89,90,92,94,95,96,98 } ;
    static const int C1_ind[] = { 1,2,3,4,6,7,8,9,12,14,16,18,22,24,26,28,31,32,33,34,36,37,38,39,42,44,46,48 };

    Matrix<double,10,10> C0; C0.setZero();
    Matrix<double,10,5> C1; C1.setZero();
    for (int i = 0; i < 66; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
    for (int i = 0; i < 28; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); }

    Matrix<double,10,5> C12 = C0.partialPivLu().solve(C1);

    // Setup action matrix
    Matrix<double,10, 5> RR;
    RR << -C12.bottomRows(5), Matrix<double,5,5>::Identity(5, 5);

    static const int AM_ind[] = { 0,1,2,3,4 };
    Matrix<double, 5, 5> AM;
    for (int i = 0; i < 5; i++) {
        AM.row(i) = RR.row(AM_ind[i]);
    }

    Matrix<std::complex<double>, 4, 5> sols;
    sols.setZero();

    // Solve eigenvalue problem
    EigenSolver<Matrix<double, 5, 5> > es(AM);
    ArrayXcd D = es.eigenvalues();
    ArrayXXcd V = es.eigenvectors();
    V = (V / V.row(4).array().replicate(5, 1)).eval();

    sols.row(0) = V.row(1).array();
    sols.row(1) = V.row(2).array();
    sols.row(2) = D.transpose().array();
    sols.row(3) = V.row(3).array();

    return sols;
}
// Action =  z
// Quotient ring basis (V) = y*w,x,y,w,1,
// Available monomials (RR*V) = y*z*w,x*z,y*z,z*w,z,y*w,x,y,w,1,
