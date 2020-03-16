#include <Eigen/Dense>
#include <complex>

using namespace Eigen;

MatrixXcd solver_floor_Hf(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(35);
    coeffs[0] = 1;
    coeffs[1] = -1;
    coeffs[2] = d[19];
    coeffs[3] = -d[0]*d[19];
    coeffs[4] = -d[1]*d[20];
    coeffs[5] = d[1]*d[19];
    coeffs[6] = d[2]*d[13] + d[3]*d[16];
    coeffs[7] = -d[0]*d[2]*d[13] - d[0]*d[3]*d[16];
    coeffs[8] = -d[1]*d[2]*d[14] - d[1]*d[3]*d[17];
    coeffs[9] = d[1]*d[2]*d[13] + d[1]*d[3]*d[16];
    coeffs[10] = -d[19];
    coeffs[11] = -d[1]*d[19];
    coeffs[12] = d[1]*d[18];
    coeffs[13] = -d[2]*d[13] - d[3]*d[16];
    coeffs[14] = -d[1]*d[2]*d[13] - d[1]*d[3]*d[16];
    coeffs[15] = d[1]*d[2]*d[12] + d[1]*d[3]*d[15];
    coeffs[16] = -d[4]*d[19];
    coeffs[17] = -d[5]*d[20];
    coeffs[18] = d[5]*d[19];
    coeffs[19] = d[6]*d[13] + d[7]*d[16];
    coeffs[20] = -d[4]*d[6]*d[13] - d[4]*d[7]*d[16];
    coeffs[21] = -d[5]*d[6]*d[14] - d[5]*d[7]*d[17];
    coeffs[22] = d[5]*d[6]*d[13] + d[5]*d[7]*d[16];
    coeffs[23] = -d[5]*d[19];
    coeffs[24] = d[5]*d[18];
    coeffs[25] = -d[6]*d[13] - d[7]*d[16];
    coeffs[26] = -d[5]*d[6]*d[13] - d[5]*d[7]*d[16];
    coeffs[27] = d[5]*d[6]*d[12] + d[5]*d[7]*d[15];
    coeffs[28] = -d[8]*d[19];
    coeffs[29] = -d[9]*d[20];
    coeffs[30] = d[9]*d[19];
    coeffs[31] = d[10]*d[13] + d[11]*d[16];
    coeffs[32] = -d[8]*d[10]*d[13] - d[8]*d[11]*d[16];
    coeffs[33] = -d[9]*d[10]*d[14] - d[9]*d[11]*d[17];
    coeffs[34] = d[9]*d[10]*d[13] + d[9]*d[11]*d[16];

	// Setup elimination template
	static const int coeffs0_ind[] = { 2,3,2,16,2,3,10,16,10,28,11,23,4,12,17,24,29,14,26,11,23,14,26,6,7,19,20,3,16,2,2,2,31,7,13,20,25,10,10,16,3,28,32,8,15,21,27,12,24,17,4,29,33,9,22,18,5,30,34 };
	static const int coeffs1_ind[] = { 7,20,19,6,31,13,25,20,7,32,15,27,21,8,33,22,9,34 };
	static const int C0_ind[] = { 0,1,2,3,9,10,11,12,13,19,21,23,30,31,32,33,39,41,43,44,45,54,55,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,92,96,97,98,99 } ;
	static const int C1_ind[] = { 4,5,6,7,8,14,15,16,17,18,24,25,26,27,28,36,37,38 };

	Matrix<double,10,10> C0; C0.setZero();
	Matrix<double,10,4> C1; C1.setZero();
	for (int i = 0; i < 59; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 18; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); }

	Matrix<double,10,4> C12 = C0.partialPivLu().solve(C1);

	// Setup action matrix
	Matrix<double,8, 4> RR;
	RR << -C12.bottomRows(4), Matrix<double,4,4>::Identity(4, 4);

	static const int AM_ind[] = { 0,1,2,3 };
	Matrix<double, 4, 4> AM;
	for (int i = 0; i < 4; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 6, 4> sols;
	sols.setZero();

	// Solve eigenvalue problem
	EigenSolver<Matrix<double, 4, 4> > es(AM);
	ArrayXcd D = es.eigenvalues();
	ArrayXXcd V = es.eigenvectors();

    // Normalize eigenvectors
    ArrayXcd normalization = (V.row(0).array().square() + V.row(1).array().square()).sqrt();

    for (int i = 0; i < 4; i++) {
        V.col(i) /= normalization(i);
    }

    sols.row(0) = V.row(0).array();
    sols.row(1) = V.row(1).array();
    sols.row(3) = V.row(2).array();
    sols.row(4) = V.row(3).array();
    sols.row(5) = D.transpose().array();

    // Extract x3 linearly from the third or fifth equations
    // It can be written as a*x_3 + b = 0, or, c*x_3 + d = 0, where
    double x7, x8, x9, x10, x11, x12, x13, x14, x19, x20, x22, x23, x25, x26;
    x7 = d[0]; x8 = d[1]; x9 = d[2]; x10 = d[3]; x11 = d[4]; x12 = d[5]; x13 = d[6]; x14 = d[7];
    x19 = d[12]; x20 = d[13]; x22 = d[15]; x23 = d[16]; x25 = d[18]; x26 = d[19];

    std::complex<double> a, b, c, dd, x1, x2, x3, x4, x6;
    for (int j=0; j < 4; j++) {
        x1 = sols(0,j);
        x2 = sols(1,j);
        x4 = sols(3,j);
        x6 = sols(5,j);
        a = -x6*x8*x26 - x8*x9*x20 - x8*x10*x23;
        c = -x6*x12*x26 - x12*x13*x20 - x12*x14*x23;
        if (std::norm(a) > std::norm(c)) {
            b = -x1*x6*x7*x26 - x1*x7*x9*x20 - x1*x7*x10*x23 - x2*x6*x26 - x2*x9*x20 - x2*x10*x23 + x4*x6*x8*x25 + x4*x8*x9*x19 + x4*x8*x10*x22;
            x3 = -b / a;
        } else {
            dd = -x1*x6*x11*x26 - x1*x11*x13*x20 - x1*x11*x14*x23 - x2*x6*x26 - x2*x13*x20 - x2*x14*x23 + x4*x6*x12*x25 + x4*x12*x13*x19 + x4*x12*x14*x22;
            x3 = -dd / c;
        }
        sols(2,j) = x3;
    }

	return sols;
}
// Action =  x6
// Quotient ring basis (V) = x1,x2,x4,x5,
// Available monomials (RR*V) = x1*x6,x2*x6,x4*x6,x5*x6,x1,x2,x4,x5,
