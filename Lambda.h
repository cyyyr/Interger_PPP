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
    static void gauss_transformation(Matrix<double> &, Matrix<double> &, int, int);

    static void permutations(Matrix<double> &, Matrix<double> &, int, double, Matrix<double> &);

    /* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) -------------------------*/
    static void reduction(Matrix<double> &, Matrix<double> &, Matrix<double> &);

    /* modified lambda (mlambda) search ---------------------------------------------------
    * args   : int    m      I  number of fixed solutions
    *          double &L     I  unit lower triangular matrix (n x n)
    *          double &D     I  diagonal matrix (n x 1)
    *          double &zs    I  (n x 1)
    *          double &zn    I  (n x m)
    *          double &s     O  sum of squared residulas of fixed solutions (1 x m)
    * return : status (0:ok,other:error)
    *
    * int    n    is a number of float parameters
    *-----------------------------------------------------------------------------------*/
    [[nodiscard]] int search(const int &, const Matrix<double> &, const Matrix<double> &,
                             Matrix<double> &, Matrix<double> &, Matrix<double> &) const;

    /* lambda/mlambda integer least-square estimation ------------------------------
    * integer least-square estimation. reduction is performed by lambda,
    * and search by mlambda.
    * args   : int    n      I  number of float parameters
    *          int    m      I  number of fixed solutions
    *          double *a     I  float parameters (n x 1)
    *          double *Q     I  covariance matrix of float parameters (n x n)
    *          double *F    O  fixed solutions (n x m)
    *          double *s     O  sum of squared residulas of fixed solutions (1 x m)
    * return : status (0:ok,other:error)
    *-----------------------------------------------------------------------------*/
    int lambda(const int &, const Matrix<double> &, Matrix<double> &,
               Matrix<double> &, Matrix<double> &);

public:
    Lambda() = default;

    ~Lambda() = default;

    Matrix<int> resolveAmbiguityWithILS(Matrix<double> &floatAmbiguity, Matrix<double>& ambiguityCovarianceMatrix);


};


#endif //ILS_PPP_LAMBDA_H