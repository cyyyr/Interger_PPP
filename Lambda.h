//
// Created by cyr on 23.10.2020.
//

#ifndef ILS_PPP_LAMBDA_H
#define ILS_PPP_LAMBDA_H

#include <vector>
#include <utility> // std::swap
#include <cmath> // std::round [(floor((x) + 0.5))]

class Lambda {
private:
    const int LOOPMAX = 10000; /* maximum count of search loop */
private:
    template<typename T>
    int sgn(T val) const {
        return (T(0) < val) - (val < T(0));
    }

private:
    static int LD_factorization(const int &n, std::vector<double> Q, std::vector<double> L, std::vector<double> D);

    static void gauss_transformation(const int &n, std::vector<double> L, std::vector<double> Z, int i, int j);

    static void
    permutations(const int &n, std::vector<double> L, std::vector<double> D, int j, double del, std::vector<double> Z);

    static void reduction(const int &n, std::vector<double> L, std::vector<double> D, const std::vector<double> &Z);

    [[nodiscard]] int search(const int &n, const int &m, const std::vector<double> &L, const std::vector<double> &D,
                             const std::vector<double> &zs,
                             std::vector<double> zn, std::vector<double> s) const;

    static int lambda_reduction(const int &n, const std::vector<double> &Q, std::vector<double> Z);

    int lambda_search(const int &n, const int &m, const std::vector<double> &a, const std::vector<double> &Q,
                      std::vector<double> F, std::vector<double> s);

public:
    int lambda(const int &n, const int &m, const std::vector<double> &a, const std::vector<double> &Q,
               const std::vector<double> &F,
               std::vector<double> s);


};


#endif //ILS_PPP_LAMBDA_H
