//
// Created by cyr on 26.10.2020.
//
#pragma once
#ifndef ILS_PPP_MATRIX_H
#define ILS_PPP_MATRIX_H

#define EPS 1e-10

#include <iostream>
#include <stdexcept>
#include <vector>

template<class T>
class Matrix {
public:
    Matrix<T>(int, int);

    [[maybe_unused]] Matrix<T>(T **, int, int);

    Matrix<T>();

    ~Matrix<T>();

    Matrix<T>(const Matrix<T> &);

    Matrix<T> &operator=(const Matrix<T> &);

    inline T &operator()(int x, int y) { return p[x][y]; }

    inline T &operator()(int x, int y) const { return p[x][y]; }

    Matrix<T> &operator+=(const Matrix<T> &);

    Matrix<T> &operator-=(const Matrix<T> &);

    Matrix<T> &operator*=(const Matrix<T> &);

    Matrix<T> &operator*=(T);

    Matrix<T> &operator/=(T);

    void swapRows(int, int);

    Matrix<T> pushBackRow(const std::vector<T> &);

    Matrix<T> pushBackColumn(const std::vector<T> &);

    Matrix<T> transpose();

    static Matrix<T> createIdentity(int);

    // functions on augmented matrices
    static Matrix<T> augment(Matrix<T>, Matrix<T>);

    Matrix<T> gaussianEliminate();

    Matrix<T> rowReduceFromGaussian();

    Matrix<T> inverse();

    [[nodiscard]] inline int getRows() const;

    [[maybe_unused]] inline void setRows(int rows);

    [[nodiscard]] inline int getCols() const;

    [[maybe_unused]] inline void setCols(int cols);

    template<typename U>
    friend std::ostream &operator<<(std::ostream &, const Matrix<U> &);

    template<typename U>
    friend std::istream &operator>>(std::istream &, Matrix<U> &);

private:
    int rows_{}, cols_{};
    T **p{};

    void allocSpace();
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

#endif //ILS_PPP_MATRIX_H

/* PUBLIC MEMBER FUNCTIONS********************************/

template<class T>
Matrix<T>::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = (T) 0;
        }
    }
}

template<class T>
[[maybe_unused]] Matrix<T>::Matrix(T **a, int rows, int cols) : rows_(rows), cols_(cols) {
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

template<class T>
Matrix<T> Matrix<T>::pushBackRow(const std::vector<T> &row) {
    if (this->cols_ == 1 && this->rows_ == 1) {
        cols_ = static_cast<int>(row.size());
        allocSpace();
        for (int i = 0; i < this->cols_; ++i) {
            this->p[rows_ - 1][i] = std::move(row[i]);
        }
    } else if (this->cols_ == static_cast<int>(row.size())) {
        ++rows_;
        T **pBuf = new T *[rows_ - 1];
        for (int i = 0; i < rows_ - 1; ++i) {
            pBuf[i] = new T[cols_];
        }
        pBuf = p;
        allocSpace();
        for (int j = 0; j < rows_ - 1; ++j) {
            for (int i = 0; i < this->cols_; ++i) {
                this->p[j][i] = std::move(pBuf[j][i]);
            }
        }
        for (int i = 0; i < this->cols_; ++i) {
            this->p[rows_ - 1][i] = std::move(row[i]);
        }
        for (int i = 0; i < rows_ - 1; ++i) {
            delete[] pBuf[i];
        }
        delete[] pBuf;
    }
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::pushBackColumn(const std::vector<T> &col) {
    if (this->cols_ == 1 && this->rows_ == 1) {
        rows_ = static_cast<int>(col.size());
        allocSpace();
        for (int i = 0; i < this->rows_; ++i) {
            this->p[i][cols_ - 1] = std::move(col[i]);
        }
    } else if (this->rows_ == static_cast<int>(col.size())) {
        ++cols_;
        T **pBuf = new T *[rows_];
        for (int i = 0; i < rows_; ++i) {
            pBuf[i] = new T[cols_ - 1];
        }
        pBuf = p;
        allocSpace();
        for (int j = 0; j < cols_ - 1; ++j) {
            for (int i = 0; i < this->rows_; ++i) {
                this->p[i][j] = std::move(pBuf[i][j]);
            }
        }
        for (int i = 0; i < this->rows_; ++i) {
            this->p[i][cols_ - 1] = std::move(col[i]);
        }
        for (int i = 0; i < rows_; ++i) {
            delete[] pBuf[i];
        }
        delete[] pBuf;
    }
    return *this;
}

/* GETTERS AND SETTERS
 ********************************/

template<class T>
inline int Matrix<T>::getRows() const {
    return rows_;
}

template<class T>
[[maybe_unused]] inline void Matrix<T>::setRows(int rows) {
    rows_ = rows;
}

template<class T>
inline int Matrix<T>::getCols() const {
    return cols_;
}

template<class T>
[[maybe_unused]] inline void Matrix<T>::setCols(int cols) {
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


/* PRIVATE HELPER FUNCTIONS
 ********************************/

template<class T>
void Matrix<T>::allocSpace() {
    p = new T *[rows_];
    for (int i = 0; i < rows_; ++i) {
        p[i] = new T[cols_];
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