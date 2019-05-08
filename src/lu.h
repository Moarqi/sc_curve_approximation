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

/// Our LU solver
class LU_Solver
{
public:

    /// empty constructor
    LU_Solver() {}

    /// factorize matrix A=L*U
    bool factorize(const MatrixXX& _A);

    /// solve A*x=b
    void solve(const VectorX& _b, VectorX& _x);

public:

    /// store factors P*A=L*U
    MatrixXX L, U, P;
};

//=============================================================================
