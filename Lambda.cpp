//
// Created by cyr on 23.10.2020.
//

#include <utility>
#include <iostream>
#include "Lambda.h"

int Lambda::factorization(const Matrix<double> &Q, Matrix<double> &L, Matrix<double> &D) {

    int n = Q.getRows();
    Matrix<double> A(Q);

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

void Lambda::gaussTransformation(Matrix<double> &L, Matrix<double> &Z, int i, int j) {

    int n = L.getRows();
    int mu = (int) std::round(L(i, j));

    if (mu != 0) {
        for (int k = i; k < n; k++)
            L(k, j) -= (double) mu * L(k, i);
        for (int k = 0; k < n; k++)
            Z(k, j) -= (double) mu * Z(k, i);
    }
}

/* TODO: try to make Z integer Matrix according to LAMBDA algorithm. */
void Lambda::permutations(Matrix<double> &L, Matrix<double> &D, int j, double del, Matrix<double> &Z) {

    const int n = L.getRows();
    double eta = D(j, 0) / del;
    double lam = D(j + 1, 0) * L(j + 1, j) / del;
    D(j, 0) = eta * D(j + 1, 0);
    D(j + 1, 0) = del;
    for (int k = 0; k <= j - 1; k++) {
        double temp1 = L(j, k);
        double temp2 = L(j + 1, k);
        L(j, k) = -L(j + 1, j) * temp1 + temp2;
        L(j + 1, k) = eta * temp1 + lam * temp2;
    }
    L(j + 1, j) = lam;
    for (int k = j + 2; k < n; k++)
        std::swap(L(k, j), L(k, j + 1));
    for (int k = 0; k < n; k++)
        std::swap(Z(k, j), Z(k, j + 1));
}

/* TODO: try to implement modified reduction. Compare speed with profiler. */
void Lambda::reduction(Matrix<double> &L, Matrix<double> &D, Matrix<double> &Z) {

    int n = L.getRows();
    int j{n - 2}, k{n - 2};

    while (j >= 0) {
        if (j <= k) {
            for (int i = j + 1; i < n; i++) {
                gaussTransformation(L, Z, i, j);
            }
        }
        double del = D(j, 0) + L(j + 1, j) * L(j + 1, j) * D(j + 1, 0);
        if (del + 1E-6 < D(j + 1, 0)) { /* compared considering numerical error */
            permutations(L, D, j, del, Z);
            k = j;
            j = n - 2;
        } else {
            j--;
        }
    }
}

int Lambda::search(const int &m, const Matrix<double> &L,
                   const Matrix<double> &D, Matrix<double> &zs,
                   Matrix<double> &zn, Matrix<double> &s) {
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

    // LD factorization
    if (!(factorization(Q, L, D))) {
        // lambda reduction
        reduction(L, D, Z);
        z = Z.transpose() * a; /* z=Z'*a */
        // mlambda
        if (!(search(m, L, D, z, E, s))) {
            F = (Z.transpose().inverse()) * E; /* F=Z'\E */
        }
    }
    return 0;
}

int Lambda::validateSolution(const Matrix<double> &S) {
    const double ambiguityFixingThreshold = 6.0; //set by user
    return (((S(0, 0) < 1e-12) ? 999.9 : S(1, 0) / S(0, 0)) < ambiguityFixingThreshold);
}

Matrix<int> Lambda::computeIntegerSolution(Matrix<double> &floatAmbiguity, Matrix<double> &ambiguityCovarianceMatrix) {
    if (floatAmbiguity.getRows() != ambiguityCovarianceMatrix.getRows() ||
        floatAmbiguity.getRows() != ambiguityCovarianceMatrix.getCols()) {
        std::cerr << "The dimensions of float ambiguity matrix and ambiguity covariance matrix does not match.\n";
    } else {
        const int m = 2; // dimensions of result matrices
        Matrix<double> F(floatAmbiguity.getRows(), m);
        Matrix<double> S(m, 1);
        if (!lambda(m, floatAmbiguity, ambiguityCovarianceMatrix, F, S)) {
            Matrix<int> ambInt(floatAmbiguity.getRows(), floatAmbiguity.getCols());
            for (int i = 0; i < floatAmbiguity.getRows(); i++) {
                ambInt(i, 0) = static_cast<int>(std::round(F(i, 0)));
            }
            if (validateSolution(S))
                return ambInt;
            else {
                std::cerr << "Integer ambiguity solution validation error\n";
                return Matrix<int>();
            }
        }
    }
    return Matrix<int>();
}
