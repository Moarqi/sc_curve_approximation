//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================
#pragma once
//=============================================================================

#include <iostream>
#include <Eigen/Dense>
#include "types.h"

//== CLASS DEFINITION =========================================================

/// Our Cholesky solver
class CholeskySolver
{
public:

    /// empty constructor
    CholeskySolver() {}

    /// factorize matrix A=L*L^T
    bool factorize(const MatrixXX& _A);

    /// solve A*x=b
    void solve(const VectorX& _b, VectorX& _x);

public:

    /// matrix factor L
    MatrixXX L;
};

//=============================================================================
