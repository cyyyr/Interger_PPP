//
// Created by cyr on 26.10.2020.
//

#ifndef ILS_PPP_MATRIX_H
#define ILS_PPP_MATRIX_H

#include <iostream>

class Matrix {
public:
    Matrix(int, int);

    Matrix(double **, int, int);

    Matrix();

    ~Matrix();

    Matrix(const Matrix &);

    Matrix &operator=(const Matrix &);

    inline double &operator()(int x, int y) { return p[x][y]; }

    Matrix &operator+=(const Matrix &);

    Matrix &operator-=(const Matrix &);

    Matrix &operator*=(const Matrix &);

    Matrix &operator*=(double);

    Matrix &operator/=(double);

    Matrix operator^(int);

    friend std::ostream &operator<<(std::ostream &, const Matrix &);

    friend std::istream &operator>>(std::istream &, Matrix &);

    void swapRows(int, int);

    Matrix transpose();

    static Matrix createIdentity(int);

    static Matrix solve(Matrix, Matrix);

    // functions on vectors
    static double dotProduct(Matrix, Matrix);

    // functions on augmented matrices
    static Matrix augment(Matrix, Matrix);

    Matrix gaussianEliminate();

    Matrix rowReduceFromGaussian();

    Matrix inverse();

private:
    int rows_, cols_;
    double **p{};

    void alloc_space();

    Matrix expHelper(const Matrix &, int);
};

Matrix operator+(const Matrix &, const Matrix &);

Matrix operator-(const Matrix &, const Matrix &);

Matrix operator*(const Matrix &, const Matrix &);

Matrix operator*(const Matrix &, double);

Matrix operator*(double, const Matrix &);

Matrix operator/(const Matrix &, double);

#endif //ILS_PPP_MATRIX_H
