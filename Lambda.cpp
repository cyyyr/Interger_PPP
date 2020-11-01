//
// Created by cyr on 23.10.2020.
//

#include <utility>
#include <iostream>
#include "Lambda.h"
#include "Matrix.h"

/* LD factorization (Q=L'*diag(D)*L) ------------------------------------------------------------------------*/
int Lambda::LD_factorization(const int &n, std::vector<double> Q, std::vector<double> L, std::vector<double> D) {
    int i, j, k;
    int info = 0;
    double a;

    std::vector<double> A(std::move(Q));

    for (i = n - 1; i >= 0; i--) {
        if ((D[i] = A[i + i * n]) <= 0.0) {
            info = -1;
            break;
        }
        a = sqrt(D[i]);
        for (j = 0; j <= i; j++) {
            L[i + j * n] = A[i + j * n] / a;
        }
        for (j = 0; j <= i - 1; j++) {
            for (k = 0; k <= j; k++) {
                A[j + k * n] -= L[i + k * n] * L[i + j * n];
            }
        }
        for (j = 0; j <= i; j++) {
            L[i + j * n] /= L[i + i * n];
        }
    }
    if (!info) {
        std::cerr << "LD factorization error: " << __FILE__ << ": "
                  << __LINE__ << '\n';
    }
    return info;
}


void Lambda::gauss_transformation(const int &n, std::vector<double> L, std::vector<double> Z, int i, int j) {
    int k;
    int mu;

    mu = static_cast<int>(std::round(L[i + j * n]));
    for (k = i; k < n; k++)
        L[k + n * j] -= mu * L[k + i * n];
    for (k = 0; k < n; k++)
        Z[k + n * j] -= mu * Z[k + i * n];
}

void
Lambda::permutations(const int &n, std::vector<double> L, std::vector<double> D, int j, double del,
                     std::vector<double> Z) {
    int k;
    double eta;
    double lam;
    double a0;
    double a1;

    eta = D[j] / del;
    lam = D[j + 1] * L[j + 1 + j * n] / del;
    D[j] = eta * D[j + 1];
    D[j + 1] = del;
    for (k = 0; k <= j - 1; k++) {
        a0 = L[j + k * n];
        a1 = L[j + 1 + k * n];
        L[j + k * n] = -L[j + 1 + j * n] * a0 + a1;
        L[j + 1 + k * n] = eta * a0 + lam * a1;
    }
    L[j + 1 + j * n] = lam;
    for (k = j + 2; k < n; k++)
        std::swap(L[k + j * n], L[k + (j + 1) * n]);
    for (k = 0; k < n; k++)
        std::swap(Z[k + j * n], Z[k + (j + 1) * n]);
}

/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
void Lambda::reduction(const int &n, std::vector<double> L, std::vector<double> D, const std::vector<double> &Z) {
    int i, j = n - 2, k = n - 2;
    double del;

    while (j >= 0) {
        if (j <= k) {
            for (i = j + 1; i < n; i++) {
                gauss_transformation(n, L, Z, i, j);
            }
        }
        del = D[j] + L[j + 1 + j * n] * L[j + 1 + j * n] * D[j + 1];
        if (del + 1E-6 < D[j + 1]) { /* compared considering numerical error */
            permutations(n, L, D, j, del, Z);
            k = j;
            j = n - 2;
        } else {
            j--;
        }
    }
}

/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
int Lambda::search(const int &n, const int &m, const std::vector<double> &L, const std::vector<double> &D,
                   const std::vector<double> &zs, std::vector<double> zn, std::vector<double> s) const {
    int i, j, k, c;
    int nn = 0, imax = 0;
    double newdist, maxdist = 1E99;
    double y;
    std::vector<double> S(n * n, 0.0);
    std::vector<double> dist(n);
    std::vector<double> zb(n);
    std::vector<double> z(n);
    std::vector<double> step(n);

    k = n - 1;
    dist[k] = 0.0;
    zb[k] = zs[k];
    z[k] = std::round(zb[k]);
    y = zb[k] - z[k];
    step[k] = static_cast<double>(sgn(y));
    for (c = 0; c < LOOPMAX; c++) {
        newdist = dist[k] + y * y / D[k];
        if (newdist < maxdist) {
            if (k != 0) {
                dist[--k] = newdist;
                for (i = 0; i <= k; i++) {
                    S[k + i * n] = S[k + 1 + i * n] + (z[k + 1] - zb[k + 1]) * L[k + 1 + i * n];
                }
                zb[k] = zs[k] + S[k + k * n];
                z[k] = std::round(zb[k]);
                y = zb[k] - z[k];
                step[k] = static_cast<double>(sgn(y));
            } else {
                if (nn < m) {
                    if (nn == 0 || newdist > s[imax]) {
                        imax = nn;
                    }
                    for (i = 0; i < n; i++) {
                        zn[i + nn * n] = z[i];
                    }
                    s[nn++] = newdist;
                } else {
                    if (newdist < s[imax]) {
                        for (i = 0; i < n; i++) {
                            zn[i + imax * n] = z[i];
                        }
                        s[imax] = newdist;
                        for (i = imax = 0; i < m; i++) {
                            if (s[imax] < s[i]) {
                                imax = i;
                            }
                        }
                    }
                    maxdist = s[imax];
                }
                z[0] += step[0];
                y = zb[0] - z[0];
                step[0] = -step[0] - static_cast<double>(sgn(step[0]));
            }
        } else {
            if (k == n - 1) {
                break;
            }

            k++;
            z[k] += step[k];
            y = zb[k] - z[k];
            step[k] = -step[k] - static_cast<double>(sgn(step[k]));
        }
    }
    for (i = 0; i < m - 1; i++) { /* sort by s */
        for (j = i + 1; j < m; j++) {
            if (s[i] < s[j]) {
                continue;
            }
            std::swap(s[i], s[j]);
            for (k = 0; k < n; k++) {
                std::swap(zn[k + i * n], zn[k + j * n]);
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
int Lambda::lambda_reduction(const int &n, const std::vector<double> &Q, std::vector<double> Z) {
    std::vector<double> L(n * n, 0.0);
    std::vector<double> D(n);
    int i;
    int j;
    int info;

    if (n <= 0) {
        return -1;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Z[i + j * n] = i == j ? 1.0 : 0.0;
        }
    }
    /* LD factorization */
    if ((info = LD_factorization(n, Q, L, D))) {
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
int Lambda::lambda_search(const int &n, const int &m, const std::vector<double> &a, const std::vector<double> &Q,
                          std::vector<double> F,
                          std::vector<double> s) {
    std::vector<double> L(n * n, 0.0);
    std::vector<double> D(n);
    int info;

    if (n <= 0 || m <= 0) {
        return -1;
    }

    /* LD factorization */
    if ((info = LD_factorization(n, Q, L, D))) {
        return info;
    }
    /* mlambda search */
    info = search(n, m, L, D, a, std::move(F), std::move(s));

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
int Lambda::lambda(const int &n, const int &m, const std::vector<double> &a, const std::vector<double> &Q,
                   const std::vector<double> &F,
                   std::vector<double> s) {
    int info;
    std::vector<double> L(n * n, 0.0);
    std::vector<double> D(n);
    std::vector<double> Z(n);
    std::vector<double> z(n);
    std::vector<double> E(n * m);

    if (n <= 0 || m <= 0) {
        return -1;
    }

    for (int i = 0; i < n; i++) {
        Z[i + i * n] = 1.0;
    }

    /* LD factorization */
    if (!(info = LD_factorization(n, Q, L, D))) {
        /* lambda reduction */
        reduction(n, L, D, Z);
        //matmul("TN", n, 1, n, 1.0, Z, a, 0.0, z); /* z=Z'*a */

        /* mlambda search */
        if (!(info = search(n, m, L, D, z, E, std::move(s)))) {
            //info = solve("T", Z, E, n, m, F); /* F=Z'\E */
        }
    }
    return info;
}
