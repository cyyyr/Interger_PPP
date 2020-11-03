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

    static void permutations(const int &, Matrix<double>, Matrix<double>, int, double, Matrix<double>);

    static void reduction(const int &, const Matrix<double> &, const Matrix<double> &, const Matrix<double> &);

    [[nodiscard]] int search(const int &, const int &, const Matrix<double> &, const Matrix<double> &,
                             const Matrix<double> &, const Matrix<double> &, const Matrix<double> &) const;

    static int lambda_reduction(const int &, const Matrix<double> &, const Matrix<double> &);

    int lambda_search(const int &, const int &, const Matrix<double> &, const Matrix<double> &,
                      const Matrix<double> &, const Matrix<double> &);

public:
    int lambda(const int &, const int &, const Matrix<double> &, const Matrix<double> &,
               Matrix<double>, const Matrix<double> &);


};


#endif //ILS_PPP_LAMBDA_H