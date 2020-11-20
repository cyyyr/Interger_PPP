//
// Created by cyr on 23.10.2020.
//

#include <utility>
#include <iostream>
#include "Lambda.h"
/*
 * TODO: 0) расположить методы в правильном порядке (см. статью)
 *       1) убрать лишнее, добавить недостающее
 *       2) добиться правильной работы программы
 *       3) написать петрову
 *       4) убрать излишек в matrix.h/cpp
 *       5) ???
*/


/* L'DL factorization (Q=L'*diag(D)*L) ------------------------------------------------------------------------*/
int Lambda::factorization(const Matrix<double> &Q, Matrix<double> &L, Matrix<double> &D) {

    int n = Q.getRows();
    Matrix<double> A(Q.getRows(), Q.getCols());
    A = Q;
//    L = Matrix<double>::Zero(n, n);
//    D = Matrix<double>::Zero(n, 1);


    for (int i = n - 1; i >= 0; i--) {
        D(i, 0) = A(i, i);
        if (D(i, 0) <= 0.0) {
            std::cerr << "LD factorization error: " << __FILE__ << ": "
                      << __LINE__ << '\n';
            return -1;
        }
        double temp = sqrt(D(i, 0));
        for (int j = 0; j <= i; j++) {
            L(i, j) = A(i, j) / temp;
        }
        for (int j = 0; j <= i - 1; j++) {
            for (int k = 0; k <= j; k++) {
                A(j, k) -= L(i, k) * L(i, j);
            }
        }
        for (int j = 0; j <= i; j++) {
            L(i, j) /= L(i, i);
        }
    }
    return 0;
}

/*Integer Gauss transformations (Z_i_j = I - mu*e_i*e'_j ; where mu is an integer ---------------*/
void Lambda::gauss_transformation(Matrix<double> &L, Matrix<double> &Z, int i, int j) {

    int n = L.getRows();
    int mu = (int) std::round(L(i, j));

    for (int k = i; k < n; k++)
        L(k, j) -= (double) mu * L(k, i);
    for (int k = 0; k < n; k++)
        Z(k, j) -= (double) mu * Z(k, i);
}


void Lambda::permutations(Matrix<double> &L, Matrix<double> &D, int j, double del, Matrix<double> &Z) {

    const int n = L.getRows();
    double eta = D(j, 0) / del;
    double lam = D(j + 1, 0) * L(j + 1, j) / del;
    D(j, 0) = eta * D(j + 1, 0);
    D(j + 1, 0) = del;
    for (int k = 0; k <= j - 1; k++) {
        L(j, k) = -L(j + 1, j) * L(j, k) + L(j + 1, k);
        L(j + 1, k) = eta * L(j, k) + lam * L(j + 1, k);
    }
    L(j + 1, j) = lam;
    for (int k = j + 2; k < n; k++)
        std::swap(L(k, j), L(k, j + 1));
    for (int k = 0; k < n; k++)
        std::swap(Z(k, j), Z(k, j + 1));
}

/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
void Lambda::reduction(Matrix<double> &L, Matrix<double> &D, Matrix<double> &Z) {

    int n = L.getRows();
    int j{n - 2}, k{n - 2};

    while (j >= 0) {
        if (j <= k) {
            for (int i = j + 1; i < n; i++) {
                gauss_transformation(L, Z, i, j);
            }
        }
        double del = D(j, 0) + L(j + 1, j) * L(j + 1, j) * D(j + 1, 0);
        std::cout << "j = " << j << ", k = " << k << ", del = " << del << '\n';
        if (del + 1E-6 < D(j + 1, 0)) { /* compared considering numerical error */
            permutations(L, D, j, del, Z);
            k = j;
            j = n - 2;
        } else {
            j--;
        }
    }
}

/* modified lambda (mlambda) search --------------------------------------------------------------
 * TODO: CHECK UNEXPECTED INPUT
 * args   : int    m      I  number of fixed solutions
 *          double &L     I  unit lower triangular matrix (n x n)
 *          double &D     I  diagonal matrix (n x 1)
 *          double &zs    I  (n x 1)
 *          double &zn    I  (n x m)
 *          double &s     O  sum of squared residulas of fixed solutions (1 x m)
 * return : status (0:ok,other:error)
 *
 * int    n    is a number of float parameters
 *-----------------------------------------------------------------------------*/
int Lambda::search(const int &m, const Matrix<double> &L,
                   const Matrix<double> &D, Matrix<double> &zs,
                   Matrix<double> &zn, Matrix<double> &s) const {
    const int maxLoopCount = 10000;
    int n = L.getRows();
    int nn{0}, imax{0}, l;
    double newDist, maxDist{1E99};
    Matrix<double> S(n, n);
    Matrix<double> dist(n, 1);
    Matrix<double> zb(n, 1);
    Matrix<double> z(n, 1);
    Matrix<double> step(n, 1);

    int k = n - 1;
    dist(k, 0) = 0.0;
    zb(k, 0) = zs(k, 0);
    z(k, 0) = std::round(zb(k, 0));
    double y = zb(k, 0) - z(k, 0);
    step(k, 0) = static_cast<double>(sign(y));
    for (l = 0; l < maxLoopCount; l++) {
        newDist = dist(k, 0) + y * y / D(k, 0);
        if (newDist < maxDist) {
            if (k != 0) {
                dist(--k, 0) = newDist;
                for (int i = 0; i <= k; i++) {
                    S(k, i) = S(k + 1, i) + (z(k + 1, 0) - zb(k + 1, 0)) * L(k + 1, 0);
                }
                zb(k, 0) = zs(k, 0) + S(k, k);
                z(k, 0) = std::round(zb(k, 0));
                y = zb(k, 0) - z(k, 0);
                step(k, 0) = static_cast<double>(sign(y));
            } else {
                if (nn < m) {
                    if (nn == 0 || newDist > s(imax, 0)) {
                        imax = nn;
                    }
                    for (int i = 0; i < n; i++) {
                        zn(i, nn) = z(i, 0);
                    }
                    s(nn++, 0) = newDist;
                } else {
                    if (newDist < s(imax, 0)) {
                        for (int i = 0; i < n; i++) {
                            zn(i, imax) = z(i, 0);
                        }
                        s(imax, 0) = newDist;
                        for (int i = imax = 0; i < m; i++) {
                            if (s(imax, 0) < s(i, 0)) {
                                imax = i;
                            }
                        }
                    }
                    maxDist = s(imax, 0);
                }
                z(0, 0) += step(0, 0);
                y = zb(0, 0) - z(0, 0);
                step(0, 0) = -step(0, 0) - static_cast<double>(sign(step(0, 0)));
            }
        } else {
            if (k == n - 1) {
                break;
            }

            k++;
            z(k, 0) += step(k, 0);
            y = zb(k, 0) - z(k, 0);
            step(k, 0) = -step(k, 0) - static_cast<double>(sign(step(k, 0)));
        }
    }
    for (int i = 0; i < m - 1; i++) { /* sort by s */
        for (int j = i + 1; j < m; j++) {
            if (s(i, 0) < s(j, 0)) {
                continue;
            }
            std::swap(s(i, 0), s(j, 0));
            for (k = 0; k < n; k++) {
                std::swap(zn(k, i), zn(k, j));
            }
        }
    }

    if (l >= maxLoopCount) {
        std::cerr << "Search loop count overflow: " << __FILE__ << ": " << __LINE__ << '\n';
        return -1;
    }
    return 0;
}

/* lambda/mlambda integer least-square estimation ------------------------------
 * integer least-square estimation. reduction is performed by lambda (ref.[1]),
 * and search by mlambda (ref.[2]).
 * args   : int    n      I  number of float parameters
 *          int    m      I  number of fixed solutions
 *          double *a     I  float parameters (n x 1)
 *          double *Q     I  covariance matrix of float parameters (n x n)
 *          integer *F    O  fixed solutions (n x m)
 *          double *s     O  sum of squared residulas of fixed solutions (1 x m)
 * return : status (0:ok,other:error)
 * notes  : matrix stored by column-major order (fortran convention)
 *-----------------------------------------------------------------------------*/
int Lambda::lambda(const int &m, const Matrix<double> &a, Matrix<double> &Q,
                   Matrix<double> &F, Matrix<double> &s) {
    if ((a.getRows() != Q.getRows()) || (Q.getRows() != Q.getCols())) return -1;
    if (m < 1) return -1;

    const int n = a.getRows();
    Matrix<double> L(n, n);
    Matrix<double> D(n, 1), z(n, 1);
    Matrix<double> E(n, m);

    Matrix<double> Z(n, n);
    Z = Matrix<double>::createIdentity(n);

    for (int i = 0; i < n; i++) {
        Z(i, i) = 1.0;
    }

    /* LD factorization */
    if (!(factorization(Q, L, D))) {
        /* lambda reduction */
        reduction(L, D, Z);
        z = Z.transpose() * a; /* z=Z'*a */
        /* mlambda search */
        if (!(search(m, L, D, z, E, s))) {
            F = (Z.transpose().inverse()) * E; /* F=Z'\E */
        }
    }
    return 0;
}

Matrix<int> Lambda::resolveIntegerAmbiguity(Matrix<double> &ambFloat, Matrix<double> &ambCov) {
    // Check input
    if (ambFloat.getRows() != ambCov.getRows() || ambFloat.getRows() != ambCov.getCols()) {
        std::cerr << "The dimension of input does not match.\n";
        std::cerr << "Cannot perform Ambiguity Resolution!\n";
    } else {
        int m = 2;
        Matrix<double> F(ambFloat.getRows(), m);
        Matrix<double> S(m, 1);
        if (!lambda(m, ambFloat, ambCov, F, S)) {
            Matrix<int> ambFixed(ambFloat.getRows(), ambFloat.getCols());
            for (int i = 0; i < ambFloat.getRows(); i++) {
                ambFixed(i, 0) = F(i, 0);
            }
            //squaredRatio = (S(0) < 1e-12) ? 9999.9 : S(1) / S(0);
            return ambFixed;
        }
    }
    return Matrix<int>();
}
