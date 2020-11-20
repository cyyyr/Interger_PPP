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
    static int factorization(const Matrix<double> &, Matrix<double> &, Matrix<double> &);

    static void gauss_transformation(Matrix<double> &, Matrix<double> &, int, int);

    static void permutations(Matrix<double> &, Matrix<double> &, int, double, Matrix<double> &);

    static void reduction(Matrix<double> &, Matrix<double> &, Matrix<double> &);

    [[nodiscard]] int search(const int &, const Matrix<double> &, const Matrix<double> &,
                             Matrix<double> &, Matrix<double> &, Matrix<double> &) const;

    int lambda(const int &, const Matrix<double> &, Matrix<double> &,
               Matrix<double> &, Matrix<double> &);

public:
    Lambda() = default;

    ~Lambda() = default;

    Matrix<int> resolveIntegerAmbiguity(Matrix<double> &ambFloat, Matrix<double>& ambCov);


};


#endif //ILS_PPP_LAMBDA_H