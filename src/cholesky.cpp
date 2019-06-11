//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#include "cholesky.h"
#include <Eigen/Dense>

//== CLASS DEFINITION =========================================================

bool CholeskySolver::factorize(const MatrixXX& A)
{
    assert(A.rows() == A.cols());

    int i, j, k;
    const int m = A.rows();

    // initialize L
    L.resize(m, m);
    L.setZero();

/**
 * \todo Compute the Cholesky factorization, i.e., compute the matrix `L`.
 * The matrix `A` is given as a parameter, matrix `L` is a member variable.
 * The matrices can be accessed by `A(i, j)` and `L(i, j)`.
 * The function should return `true` on success and `false` on failure.
 */
    L = A.triangularView<Eigen::UpLoType::Lower>();

    for (k = 0; k < m; k++) {
        float L_kk = L(k, k);
        float one_div_by_sqrt_L_kk = 1.0 / sqrt(L_kk);
        for (j = k + 1; j < m; j++) {
            float L_jk = L(j, k);
            for (i = j; i < m; i++) {
                L(i, j) = L(i, j) - L(i, k) * L_jk / L_kk;
            }
        }
        for (i = k; i < m; i++) {
            L(i, k) = L(i, k) * one_div_by_sqrt_L_kk;
        }
    }

    // check error of factorization
    std::cout << "  error(A=L*L^T) : " << (A-L*L.transpose()).norm()
              << std::endl;

    return true;
}

//-----------------------------------------------------------------------------

void CholeskySolver::solve(const VectorX& _b, VectorX& _x)
{
    const int m = L.rows();
    assert((int)_b.size() == m);
    assert((int)_x.size() == m);

    int i, j;
    Scalar sum;

    // solve L*y=b by forward substitution
    VectorX y(m);
    for (i = 0; i < m; ++i)
    {
        sum = 0.0;
        for (j = 0; j < i; ++j) sum += L(i, j) * y(j);

        y(i) = (_b(i) - sum) / L(i, i);
    }

    // solve L^T*x=y by backward substitution
    for (i = m - 1; i >= 0; --i)
    {
        sum = 0.0;
        for (j = i + 1; j < m; ++j) sum += L(j, i) * _x(j);

        _x(i) = (y(i) - sum) / L(i, i);
    }
}

//=============================================================================
