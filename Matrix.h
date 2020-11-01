//
// Created by cyr on 26.10.2020.
//
#pragma once
#ifndef ILS_PPP_MATRIX_H
#define ILS_PPP_MATRIX_H

#define EPS 1e-10

#include <iostream>
#include <stdexcept>

template<class T>
class Matrix {
public:
    Matrix<T>(int, int);

    Matrix<T>(T **, int, int);

    Matrix<T>();

    ~Matrix<T>();

    Matrix<T>(const Matrix<T> &);

    Matrix<T> &operator=(const Matrix<T> &);

    inline T &operator()(int x, int y) { return p[x][y]; }

    Matrix<T> &operator+=(const Matrix<T> &);

    Matrix<T> &operator-=(const Matrix<T> &);

    Matrix<T> &operator*=(const Matrix<T> &);

    Matrix<T> &operator*=(T);

    Matrix<T> &operator/=(T);

    Matrix<T> operator^(int);

    void swapRows(int, int);

    Matrix<T> transpose();

    static Matrix<T> createIdentity(int);

    static Matrix<T> solve(Matrix<T>, Matrix<T>);

    // functions on vectors
    static T dotProduct(Matrix<T>, Matrix<T>);

    // functions on augmented matrices
    static Matrix<T> augment(Matrix<T>, Matrix<T>);

    Matrix<T> gaussianEliminate();

    Matrix<T> rowReduceFromGaussian();

    Matrix<T> inverse();

    int LD_factorization(const int &, Matrix<T>, Matrix<T>, Matrix<T>);

    [[nodiscard]] int getRows() const;

    void setRows(int rows);

    [[nodiscard]] int getCols() const;

    void setCols(int cols);

private:
    int rows_, cols_;
    T **p{};

    void allocSpace();

    Matrix<T> expHelper(const Matrix<T> &, int);
};

template<class T>
Matrix<T> operator+(const Matrix<T> &, const Matrix<T> &);

template<class T>
Matrix<T> operator-(const Matrix<T> &, const Matrix<T> &);

template<class T>
Matrix<T> operator*(const Matrix<T> &, const Matrix<T> &);

template<class T>
Matrix<T> operator*(const Matrix<T> &, T);

template<class T>
Matrix<T> operator*(T, const Matrix<T> &);

template<class T>
Matrix<T> operator/(const Matrix<T> &, T);

template<class T>
std::ostream &operator<<(std::ostream &, const Matrix<T> &);

template<class T>
std::istream &operator>>(std::istream &, Matrix<T> &);

#endif //ILS_PPP_MATRIX_H

/* PUBLIC MEMBER FUNCTIONS********************************/

template<class T>
Matrix<T>::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = 0;
        }
    }
}

template<class T>
Matrix<T>::Matrix(T **a, int rows, int cols) : rows_(rows), cols_(cols) {
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = a[i][j];
        }
    }
}

template<class T>
Matrix<T>::Matrix() : rows_(1), cols_(1) {
    allocSpace();
    p[0][0] = 0;
}

template<class T>
Matrix<T>::~Matrix() {
    for (int i = 0; i < rows_; ++i) {
        delete[] p[i];
    }
    delete[] p;
}

template<class T>
Matrix<T>::Matrix(const Matrix &m) : rows_(m.rows_), cols_(m.cols_) {
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
}


template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &m) {
    if (this == &m) {
        return *this;
    }

    if (rows_ != m.rows_ || cols_ != m.cols_) {
        for (int i = 0; i < rows_; ++i) {
            delete[] p[i];
        }
        delete[] p;

        rows_ = m.rows_;
        cols_ = m.cols_;
        allocSpace();
    }

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix &m) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] += m.p[i][j];
        }
    }
    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix &m) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] -= m.p[i][j];
        }
    }
    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator*=(const Matrix &m) {
    Matrix temp(rows_, m.cols_);
    for (int i = 0; i < temp.rows_; ++i) {
        for (int j = 0; j < temp.cols_; ++j) {
            for (int k = 0; k < cols_; ++k) {
                temp.p[i][j] += (p[i][k] * m.p[k][j]);
            }
        }
    }
    return (*this = temp);
}

template<class T>
Matrix<T> &Matrix<T>::operator*=(T num) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] *= num;
        }
    }
    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator/=(T num) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] /= num;
        }
    }
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator^(int num) {
    Matrix temp(*this);
    return expHelper(temp, num);
}

template<class T>
void Matrix<T>::swapRows(int r1, int r2) {
    double *temp = p[r1];
    p[r1] = p[r2];
    p[r2] = temp;
}

template<class T>
Matrix<T> Matrix<T>::transpose() {
    Matrix matrix(cols_, rows_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            matrix.p[j][i] = p[i][j];
        }
    }
    return matrix;
}

/* GETTERS AND SETTERS
 ********************************/

template<class T>
int Matrix<T>::getRows() const {
    return rows_;
}

template<class T>
void Matrix<T>::setRows(int rows) {
    rows_ = rows;
}

template<class T>
int Matrix<T>::getCols() const {
    return cols_;
}

template<class T>
void Matrix<T>::setCols(int cols) {
    cols_ = cols;
}

/* STATIC CLASS FUNCTIONS
 ********************************/

template<class T>
Matrix<T> Matrix<T>::createIdentity(int size) {
    Matrix<T> temp(size, size);
    for (int i = 0; i < temp.rows_; ++i) {
        for (int j = 0; j < temp.cols_; ++j) {
            if (i == j) {
                temp.p[i][j] = 1;
            } else {
                temp.p[i][j] = 0;
            }
        }
    }
    return temp;
}

template<class T>
Matrix<T> Matrix<T>::solve(Matrix A, Matrix b) { //A*x=b
    // Gaussian elimination
    for (int i = 0; i < A.rows_; ++i) {
        if (A.p[i][i] == 0) {
            // pivot - 0 -> throw error
            throw std::domain_error(
                    "Error: the coefficient matrix has 0 as a pivot. Please fix the input and try again.");
        }
        for (int j = i + 1; j < A.rows_; ++j) {
            for (int k = i + 1; k < A.cols_; ++k) {
                A.p[j][k] -= A.p[i][k] * (A.p[j][i] / A.p[i][i]);
                if (A.p[j][k] < EPS && A.p[j][k] > -1 * EPS)
                    A.p[j][k] = 0;
            }
            b.p[j][0] -= b.p[i][0] * (A.p[j][i] / A.p[i][i]);
            if (A.p[j][0] < EPS && A.p[j][0] > -1 * EPS)
                A.p[j][0] = 0;
            A.p[j][i] = 0;
        }
    }

    // Back substitution
    Matrix<T> x(b.rows_, 1);
    x.p[x.rows_ - 1][0] = b.p[x.rows_ - 1][0] / A.p[x.rows_ - 1][x.rows_ - 1];
    if (x.p[x.rows_ - 1][0] < EPS && x.p[x.rows_ - 1][0] > -1 * EPS)
        x.p[x.rows_ - 1][0] = 0;
    for (int i = x.rows_ - 2; i >= 0; --i) {
        T sum = 0.0;
        for (int j = i + 1; j < x.rows_; ++j) {
            sum += A.p[i][j] * x.p[j][0];
        }
        x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
        if (x.p[i][0] < EPS && x.p[i][0] > -1 * EPS)
            x.p[i][0] = 0;
    }

    return x;
}

// functions on VECTORS
template<class T>
[[maybe_unused]] T Matrix<T>::dotProduct(Matrix a, Matrix b) {
    double sum = 0;
    for (int i = 0; i < a.rows_; ++i) {
        sum += (a(i, 0) * b(i, 0));
    }
    return sum;
}

// functions on AUGMENTED matrices
template<class T>
Matrix<T> Matrix<T>::augment(Matrix A, Matrix B) {
    Matrix AB(A.rows_, A.cols_ + B.cols_);
    for (int i = 0; i < AB.rows_; ++i) {
        for (int j = 0; j < AB.cols_; ++j) {
            if (j < A.cols_)
                AB(i, j) = A(i, j);
            else
                AB(i, j) = B(i, j - B.cols_);
        }
    }
    return AB;
}

template<class T>
Matrix<T> Matrix<T>::gaussianEliminate() {
    Matrix Ab(*this);
    int rows = Ab.rows_;
    int cols = Ab.cols_;
    int Acols = cols - 1;

    int i = 0; // row tracker
    int j = 0; // column tracker

    // iterate through the rows
    while (i < rows) {
        // find a pivot for the row
        bool pivot_found = false;
        while (j < Acols && !pivot_found) {
            if (Ab(i, j) != 0) { // pivot not equal to 0
                pivot_found = true;
            } else { // check for a possible swap
                int max_row = i;
                double max_val = 0;
                for (int k = i + 1; k < rows; ++k) {
                    double cur_abs = Ab(k, j) >= 0 ? Ab(k, j) : -1 * Ab(k, j);
                    if (cur_abs > max_val) {
                        max_row = k;
                        max_val = cur_abs;
                    }
                }
                if (max_row != i) {
                    Ab.swapRows(max_row, i);
                    pivot_found = true;
                } else {
                    j++;
                }
            }
        }

        // perform elimination as normal if pivot was found
        if (pivot_found) {
            for (int t = i + 1; t < rows; ++t) {
                for (int s = j + 1; s < cols; ++s) {
                    Ab(t, s) = Ab(t, s) - Ab(i, s) * (Ab(t, j) / Ab(i, j));
                    if (Ab(t, s) < EPS && Ab(t, s) > -1 * EPS)
                        Ab(t, s) = 0;
                }
                Ab(t, j) = 0;
            }
        }

        i++;
        j++;
    }

    return Ab;
}

template<class T>
Matrix<T> Matrix<T>::rowReduceFromGaussian() {
    Matrix<T> R(*this);
    int rows = R.rows_;
    int cols = R.cols_;

    int i = rows - 1; // row tracker
    int j = cols - 2; // column tracker

    // iterate through every row
    while (i >= 0) {
        // find the pivot column
        int k = j - 1;
        while (k >= 0) {
            if (R(i, k) != 0)
                j = k;
            k--;
        }

        // zero out elements above pivots if pivot not 0
        if (R(i, j) != 0) {

            for (int t = i - 1; t >= 0; --t) {
                for (int s = 0; s < cols; ++s) {
                    if (s != j) {
                        R(t, s) = R(t, s) - R(i, s) * (R(t, j) / R(i, j));
                        if (R(t, s) < EPS && R(t, s) > -1 * EPS)
                            R(t, s) = 0;
                    }
                }
                R(t, j) = 0;
            }

            // divide row by pivot
            for (k = j + 1; k < cols; ++k) {
                R(i, k) = R(i, k) / R(i, j);
                if (R(i, k) < EPS && R(i, k) > -1 * EPS)
                    R(i, k) = 0;
            }
            R(i, j) = 1;

        }

        i--;
        j--;
    }

    return R;
}

template<class T>
Matrix<T> Matrix<T>::inverse() {
    Matrix<T> I = Matrix::createIdentity(rows_);
    Matrix<T> AI = Matrix::augment(*this, I);
    Matrix<T> U = AI.gaussianEliminate();
    Matrix<T> IAInverse = U.rowReduceFromGaussian();
    Matrix<T> AInverse(rows_, cols_);
    for (int i = 0; i < AInverse.rows_; ++i) {
        for (int j = 0; j < AInverse.cols_; ++j) {
            AInverse(i, j) = IAInverse(i, j + cols_);
        }
    }
    return AInverse;
}


/* LD factorization (Q=L'*diag(D)*L) ------------------------------------------------------------------------*/
template<class T>
int Matrix<T>::LD_factorization(const int &n, Matrix<T> Q, Matrix<T> L, Matrix<T> D) {
    int i, j, k;
    int info = 0;
    T a;

    Matrix<T> A(std::move(Q));

    for (i = n - 1; i >= 0; i--) {
        if ((D(i, 1) = A(i, i)) <= 0.0) {
            info = -1;
            break;
        }
        a = sqrt(D(i, 1));
        for (j = 0; j <= i; j++) {
            L(i, j) = A(i, j) / a;
        }
        for (j = 0; j <= i - 1; j++) {
            for (k = 0; k <= j; k++) {
                A(j, k) -= L(i, k) * L(i, j);
            }
        }
        for (j = 0; j <= i; j++) {
            L(i, j) /= L(i, i);
        }
    }
    if (!info) {
        std::cerr << "LD factorization error: " << __FILE__ << ": "
                  << __LINE__ << '\n';
    }
    return info;
}


/* PRIVATE HELPER FUNCTIONS
 ********************************/

template<class T>
void Matrix<T>::allocSpace() {
    p = new T *[rows_];
    for (int i = 0; i < rows_; ++i) {
        p[i] = new T[cols_];
    }
}

template<class T>
Matrix<T> Matrix<T>::expHelper(const Matrix &m, int num) {
    if (num == 0) {
        return createIdentity(m.rows_);
    } else if (num == 1) {
        return m;
    } else if (num % 2 == 0) {  // num is even
        return expHelper(m * m, num / 2);
    } else {                    // num is odd
        return m * expHelper(m * m, (num - 1) / 2);
    }
}

/* NON-MEMBER FUNCTIONS
 ********************************/

template<class T>
Matrix<T> operator+(const Matrix<T> &m1, const Matrix<T> &m2) {
    Matrix temp(m1);
    return (temp += m2);
}

template<class T>
Matrix<T> operator-(const Matrix<T> &m1, const Matrix<T> &m2) {
    Matrix<T> temp(m1);
    return (temp -= m2);
}

template<class T>
Matrix<T> operator*(const Matrix<T> &m1, const Matrix<T> &m2) {
    Matrix<T> temp(m1);
    return (temp *= m2);
}

template<class T>
Matrix<T> operator*(const Matrix<T> &m, const T &num) {
    Matrix<T> temp(m);
    return (temp *= num);
}

template<class T>
Matrix<T> operator*(const T &num, const Matrix<T> &m) {
    return (m * num);
}

template<class T>
Matrix<T> operator/(const Matrix<T> &m, const T &num) {
    Matrix<T> temp(m);
    return (temp /= num);
}

template<class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &m) {
    for (int i = 0; i < m.rows_; ++i) {
        os << m.p[i][0];
        for (int j = 1; j < m.cols_; ++j) {
            os << " " << m.p[i][j];
        }
        os << std::endl;
    }
    return os;
}

template<class T>
std::istream &operator>>(std::istream &is, Matrix<T> &m) {
    for (int i = 0; i < m.rows_; ++i) {
        for (int j = 0; j < m.cols_; ++j) {
            is >> m.p[i][j];
        }
    }
    return is;
}