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
    const int LOOPMAX = 10000; /* maximum count of search loop */
private:
    template<typename T>
    int sgn(T val) const {
        return (T(0) < val) - (val < T(0));
    }

private:
    static void gauss_transformation(const int &n, Matrix<double> L, Matrix<double> Z, int i, int j);

    static void
    permutations(const int &n, Matrix<double> L, Matrix<double> D, int j, double del,
                 Matrix<double> Z);

    static void reduction(const int &n, Matrix<double> L, Matrix<double> D, const Matrix<double> &Z);

    [[nodiscard]] int search(const int &n, const int &m, const Matrix<double> &L, const Matrix<double> &D,
                             const Matrix<double> &zs,
                             Matrix<double> zn, Matrix<double> s) const;

    static int lambda_reduction(const int &n, const Matrix<double> &Q, Matrix<double> Z);

    int lambda_search(const int &n, const int &m, const Matrix<double> &a, const Matrix<double> &Q,
                      Matrix<double> F, Matrix<double> s);

public:
    int lambda(const int &n, const int &m, const Matrix<double> &a, const Matrix<double> &Q,
               Matrix<double> F, const Matrix<double> &s);


};


#endif //ILS_PPP_LAMBDA_H