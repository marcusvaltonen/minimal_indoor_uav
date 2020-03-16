#include <iostream>
#include <Eigen/Dense>
#include <chrono>

#include "get_floor_fHf.h"
#include "get_floor_Hf.h"
#include "get_floor_f1Hf2.h"

using namespace Eigen;
using namespace std;

int main()
{
    /* Timing experiments */
    int nbr_iter = 10000;
    MatrixXcd out;

    // f1Hf2
    cout << "f1Hf2" << endl;
    auto start = chrono::steady_clock::now();
    for (int i=0; i < nbr_iter; i++) {
        const VectorXd data = VectorXd::Random(45);
        out = solver_floor_f1Hf2(data);
    }

    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in microseconds : "
        << chrono::duration_cast<chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << endl;

    // fHf
    cout << "fHf" << endl;
    start = chrono::steady_clock::now();
    for (int i=0; i < nbr_iter; i++) {
        const VectorXd data = VectorXd::Random(3 * 2 * 2 + 9 * 2);
        out = solver_floor_fHf(data);
    }

    end = chrono::steady_clock::now();
    cout << "Elapsed time in microseconds : "
        << chrono::duration_cast<chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << endl;

    // Hf
    cout << "Hf" << endl;
    start = chrono::steady_clock::now();
    for (int i=0; i < nbr_iter; i++) {
        const VectorXd data = VectorXd::Random(3 * 2 * 2 + 9);
        out = solver_floor_Hf(data);
    }

    end = chrono::steady_clock::now();
    cout << "Elapsed time in microseconds : "
        << chrono::duration_cast<chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << endl;

    return 0;
}
