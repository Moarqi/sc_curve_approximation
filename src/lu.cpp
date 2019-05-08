//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#include "lu.h"

//== CLASS DEFINITION =========================================================

bool LU_Solver::factorize(const MatrixXX& A)
{
    assert(A.rows() == A.cols());

    if (A.rows() < 1)
        return true;

    int i, j, k;
    const int m = A.rows();

    // initialize L, U
    L.resize(m, m);
    L = MatrixXX::Identity(m, m);
    U.resize(m, m);
    U = A;

    // main loop
    for (k = 0; k < m; k++)
    {
        std::cout << k << std::endl;
        // if the next diagonal element is too small, the matrix is singular
        if (fabs(U(k, k)) < 1e-20)
        {
            std::cerr << "LU_Solver: Factorization failed." << std::endl;
            return false;
        }
        // find pivot element of matrix
        // can i use EIGEN's implementation for that?

        // double temp_col_max = 0.0;
        // int index_col_max = 0;
        // // do it manual
        // for (i = k; i < m; i++) {
        //     if (U(k, i) > temp_col_max) {
        //         temp_col_max = U(i, k);
        //         index_col_max = i;
        //     }
        // }
        for (i = k+1; i < m; i++) {
            L(i, k) = U(i, k) / U(k, k);
            for (j = k; j < m; j++) {
                U(i, j) = U(i, j) - L(i, k) * U(k, j);
            }
        }
    }

    // check error of LU factorization
    std::cout << "  error(A=L*U) : " << (A-L*U).norm() << std::endl;

    return true;
}

//-----------------------------------------------------------------------------

void LU_Solver::solve(const VectorX& _b, VectorX& _x)
{
/**
 * \todo Solve the system `A * _x = b`, using the computed factorization of the
 * matrix `A = L * U`.
 * The right hand side is `_b`, the result is to be written to `_x`.
 * - Solve `L * y = b`, then check the error `norm(L * y - b)`.
 * - Solve `U * x = y`, then check the error `norm(U * x - y)`.
 */

    assert(_b.rows() == _x.rows());

    if (_b.rows() < 1) return;

    int i, m, j;

    m = _b.rows();

    VectorX y = VectorX(m);
    std::cout << "m: " << m << std::endl;

    std::cout << "U_rows: " << U.rows() << " U_columns: " << U.cols() << std::endl;

    std::cout << "L_rows: " << L.rows() << " L_columns: " << L.cols() << std::endl;
    y(0) = _b(0) / L(0, 0);         

    for (i = 1; i < m; i++) {
        double sum_so_far = 0;
        for (j = 0; j < i; j++) { // TODO: this might easier if I would use vector multiplication
            // multipliziere i-te row von L mit y (alles andere sollte 0 sein)
            sum_so_far += y(j) * L(i, j);
        }
        y(i) = (_b(i) - sum_so_far) / L(i, i);        
    }

    _x(m-1) = y(m-1) / U(m-1, m-1);

    for (i = m - 2; i >= 0; i--) {
        double sum_so_far = 0;
        for (j = i+1; j < m; j++) { // TODO: this might easier if I would use vector multiplication            
            sum_so_far += _x(j) * U(i, j);
        }
        _x(i) = (y(i) - sum_so_far) / U(i, i);        
    }


}

//=============================================================================
