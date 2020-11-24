//
// Created by cyr on 23.10.2020.
//
#pragma once
#ifndef ILS_PPP_LAMBDA_H
#define ILS_PPP_LAMBDA_H

#include <vector>
#include <utility> // std::swap
#include <cmath> // std::round [(floor((x) + 0.5))]
#include "Matrix.h"

class Lambda {
private:
    template<typename T>
    int sign(T val) const {
        return (T(0) < val) - (val <= T(0));
    }

private:
    /* L'DL factorization (Q=L'*diag(D)*L) ------------------------------------------------*/
    static int factorization(const Matrix<double> &, Matrix<double> &, Matrix<double> &);

    /*Integer Gauss transformations (Z_i_j = I - mu*e_i*e'_j ; where mu is an integer -----*/
    static void gaussTransformation(Matrix<double> &, Matrix<double> &, int, int);

    static void permutations(Matrix<double> &, Matrix<double> &, int, double, Matrix<double> &);

    /* lambda reduction (z = Z' * a, Qz = Z' * Q * Z = L' * diag(D) * L) --------------*/
    static void reduction(Matrix<double> &, Matrix<double> &, Matrix<double> &);

    /* mlambda search -------------------------------------------------------------------
    * input  : const int &m     number of fixed solutions
    *          const double &L  unit lower triangular matrix [n x n]
    *          const double &D  diagonal matrix [n x 1]
    *          double &zs       [n x 1]
    *          double &zn       [n x m]
    * output:  double &s        sum of squared residuals of fixed solutions [m x 1]
    * return : status : error if not 0
    *
    * int    n    is a number of float parameters
    *-----------------------------------------------------------------------------------*/
    [[nodiscard]] int search(const int &, const Matrix<double> &, const Matrix<double> &,
                             Matrix<double> &, Matrix<double> &, Matrix<double> &);

    /* lambda/mlambda integer least-square estimation ------------------------------
    * integer least-square estimation. reduction is performed by lambda,
    * and search by mlambda.
    * input  : const int &m      number of fixed solutions
    *          const double &a   float parameters [n x 1]
    *          const double &Q   covariance matrix of float parameters [n x n]
    * output :       double &F   fixed solutions [n x m]
    *                double &s   sum of squared residuals of fixed solutions [m x 1]
    * return : status : error if not 0
    *-----------------------------------------------------------------------------*/
    [[nodiscard]] int lambda(const int &, const Matrix<double> &, Matrix<double> &,
                             Matrix<double> &, Matrix<double> &);

    /* validation of the fixed solution */
    static int validateSolution(const Matrix<double> &S);

public:
    Lambda() = default;

    ~Lambda() = default;


    Matrix<int> computeIntegerSolution(Matrix<double> &floatAmbiguity, Matrix<double> &ambiguityCovarianceMatrix);


};


#endif //ILS_PPP_LAMBDA_H