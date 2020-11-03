//
// Created by cyr on 23.10.2020.
//

#include <utility>
#include <iostream>
#include "Lambda.h"

void Lambda::gauss_transformation(const int &n, Matrix<double> L, Matrix<double> Z, int i, int j) {
    int k;
    int mu;

    mu = static_cast<int>(std::round(L(i, j)));
    for (k = i; k < n; k++)
        L(k, j) -= mu * L(k, i);
    for (k = 0; k < n; k++)
        Z(k, j) -= mu * Z(k, i);
}

void Lambda::permutations(const int &n, Matrix<double> L, Matrix<double> D, int j, double del,
                          Matrix<double> Z) {
    int k;
    double eta;
    double lam;
    double a0;
    double a1;

    eta = D(j, 1) / del;
    lam = D(j + 1, 1) * L(j + 1, j) / del;
    D(j, 1) = eta * D(j + 1, 1);
    D(j + 1, 1) = del;
    for (k = 0; k <= j - 1; k++) {
        a0 = L(j, k);
        a1 = L(j + 1, k);
        L(j, k) = -L(j + 1, j) * a0 + a1;
        L(j + 1, k) = eta * a0 + lam * a1;
    }
    L(j + 1, j) = lam;
    for (k = j + 2; k < n; k++)
        std::swap(L(k, j), L(k, j + 1));
    for (k = 0; k < n; k++)
        std::swap(Z(k, j), Z(k, j + 1));
}

/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
void Lambda::reduction(const int &n, const Matrix<double> &L, const Matrix<double> &D, const Matrix<double> &Z) {
    int i, j = n - 2, k = n - 2;
    double del;

    while (j >= 0) {
        if (j <= k) {
            for (i = j + 1; i < n; i++) {
                gauss_transformation(n, L, Z, i, j);
            }
        }
        del = D(j, 1) + L(j + 1, j) * L(j + 1, j) * D(j + 1, 1);
        if (del + 1E-6 < D(j + 1, 1)) { /* compared considering numerical error */
            permutations(n, L, D, j, del, Z);
            k = j;
            j = n - 2;
        } else {
            j--;
        }
    }
}

/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
int Lambda::search(const int &n, const int &m, const Matrix<double> &L, const Matrix<double> &D,
                   const Matrix<double> &zs, const Matrix<double> &zn, const Matrix<double> &s) const {
    int i, j, k, c;
    int nn = 0, imax = 0;
    double newdist, maxdist = 1E99;
    double y;
    Matrix<double> S(n, n);
    Matrix<double> dist(n, 1);
    Matrix<double> zb(n, 1);
    Matrix<double> z(n, 1);
    Matrix<double> step(n, 1);

    k = n - 1;
    dist(k, 1) = 0.0;
    zb(k, 1) = zs(k, 1);
    z(k, 1) = std::round(zb(k, 1));
    y = zb(k, 1) - z(k, 1);
    step(k, 1) = static_cast<double>(sgn(y));
    for (c = 0; c < LOOPMAX; c++) {
        newdist = dist(k, 1) + y * y / D(k, 1);
        if (newdist < maxdist) {
            if (k != 0) {
                dist(--k, 1) = newdist;
                for (i = 0; i <= k; i++) {
                    S(k, i) = S(k + 1, i) + (z(k + 1, 1) - zb(k + 1, 1)) * L(k + 1, 1);
                }
                zb(k,1) = zs(k,1) + S(k,k);
                z(k,1) = std::round(zb(k,1));
                y = zb(k,1) - z(k,1);
                step(k,1) = static_cast<double>(sgn(y));
            } else {
                if (nn < m) {
                    if (nn == 0 || newdist > s(imax,1)) {
                        imax = nn;
                    }
                    for (i = 0; i < n; i++) {
                        zn(i,nn) = z(i,1);
                    }
                    s(nn++,1) = newdist;
                } else {
                    if (newdist < s(imax,1)) {
                        for (i = 0; i < n; i++) {
                            zn(i,imax) = z(i,1);
                        }
                        s(imax,1) = newdist;
                        for (i = imax = 0; i < m; i++) {
                            if (s(imax,1) < s(i,1)) {
                                imax = i;
                            }
                        }
                    }
                    maxdist = s(imax,1);
                }
                z(0,1) += step(0,1);
                y = zb(0,1) - z(0,1);
                step(0,1) = -step(0,1) - static_cast<double>(sgn(step(0,1)));
            }
        } else {
            if (k == n - 1) {
                break;
            }

            k++;
            z(k,1) += step(k,1);
            y = zb(k,1) - z(k,1);
            step(k,1) = -step(k,1) - static_cast<double>(sgn(step(k,1)));
        }
    }
    for (i = 0; i < m - 1; i++) { /* sort by s */
        for (j = i + 1; j < m; j++) {
            if (s(i,1) < s(j,1)) {
                continue;
            }
            std::swap(s(i,1), s(j,1));
            for (k = 0; k < n; k++) {
                std::swap(zn(k,i), zn(k,j));
            }
        }
    }

    if (c >= LOOPMAX) {
        std::cerr << "Search loop count overflow: " << __FILE__ << ": " << __LINE__ << '\n';
        return -1;
    }
    return 0;
}

/* lambda reduction ------------------------------------------------------------
 * reduction by lambda (ref [1]) for integer least square
 * args   : int    n      I  number of float parameters
 *          double *Q     I  covariance matrix of float parameters (n x n)
 *          double *Z     O  lambda reduction matrix (n x n)
 * return : status (0:ok,other:error)
 *-----------------------------------------------------------------------------*/
int Lambda::lambda_reduction(const int &n, const Matrix<double> &Q, const Matrix<double> &Z) {
    Matrix<double> L(n, n);
    Matrix<double> D(n, 1);
    int i;
    int j;
    int info;

    if (n <= 0) {
        return -1;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Z(i, j) = i == j ? 1.0 : 0.0;
        }
    }
    /* LD factorization */
    if ((info = Matrix<double>::LD_factorization(n, Q, L, D))) {
        return info;
    }
    /* lambda reduction */
    reduction(n, L, D, Z);

    return 0;
}

/* mlambda search --------------------------------------------------------------
 * search by  mlambda (ref [2]) for integer least square
 * args   : int    n      I  number of float parameters
 *          int    m      I  number of fixed solutions
 *          double *a     I  float parameters (n x 1)
 *          double *Q     I  covariance matrix of float parameters (n x n)
 *          double *F     O  fixed solutions (n x m)
 *          double *s     O  sum of squared residulas of fixed solutions (1 x m)
 * return : status (0:ok,other:error)
 *-----------------------------------------------------------------------------*/
int Lambda::lambda_search(const int &n, const int &m, const Matrix<double> &a, const Matrix<double> &Q,
                          const Matrix<double> &F, const Matrix<double> &s) {
    Matrix<double> L(n, n);
    Matrix<double> D(n,1);
    int info;

    if (n <= 0 || m <= 0) {
        return -1;
    }

    /* LD factorization */
    if ((info = Matrix<double>::LD_factorization(n, Q, L, D))) {
        return info;
    }
    /* mlambda search */
    info = search(n, m, L, D, a, F, s);

    return info;
}

/* lambda/mlambda integer least-square estimation ------------------------------
 * integer least-square estimation. reduction is performed by lambda (ref.[1]),
 * and search by mlambda (ref.[2]).
 * args   : int    n      I  number of float parameters
 *          int    m      I  number of fixed solutions
 *          double *a     I  float parameters (n x 1)
 *          double *Q     I  covariance matrix of float parameters (n x n)
 *          double *F     O  fixed solutions (n x m)
 *          double *s     O  sum of squared residulas of fixed solutions (1 x m)
 * return : status (0:ok,other:error)
 * notes  : matrix stored by column-major order (fortran convention)
 *-----------------------------------------------------------------------------*/
int Lambda::lambda(const int &n, const int &m, const Matrix<double> &a, const Matrix<double> &Q,
                   Matrix<double> F, const Matrix<double> &s) {
    if (n <= 0 || m <= 0) {
        return -1;
    }
    int info;
    Matrix<double> L(n, n);
    Matrix<double> D(n, 1);
    Matrix<double> Z(n, 1);
    Matrix<double> z(n, 1);
    Matrix<double> E(n, m);

    for (int i = 0; i < n; ++i) {
        Z(i, i) = 1.0;
    }

    /* LD factorization */
    if (!(info = Matrix<double>::LD_factorization(n, Q, L, D))) {
        /* lambda reduction */
        reduction(n, L, D, Z);
        z = Z.transpose() * a; /* z=Z'*a */
        /* mlambda search */
        if (!(info = search(n, m, L, D, z, E, s))) {
            F = Matrix<double>::solve(Z.transpose(), E); /* F=Z'\E */
        }
    }
    return info;
}
