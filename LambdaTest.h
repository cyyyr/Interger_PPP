#pragma once
#ifndef ILS_PPP_LAMBDATEST_H
#define ILS_PPP_LAMBDATEST_H

#include <istream>
#include <vector>

struct LambdaTest {
LambdaTest() = default;
~LambdaTest() = default;

    static int test();

    static int testFromFile();

    //TODO: to implement this:
//    static inline void readInputFromFile(std::ifstream &F, std::vector<double>, std::vector<double>);
//    inline void LambdaTest::readInputFromFile(std::ifstream &F, std::vector<double> inputAmbiguity,
//                                              std::vector<double> inputFloatAmbCovariance) {
//        int n{0};
//        F >> n;
//        if (n != 0) {
//            for (int i = 0; i < n; ++i) {
//                inputAmbiguity.push_back(0.0);
//                F >> inputAmbiguity[i];
//            }
//            for (int i = 0; i < n * n; ++i) {
//                inputFloatAmbCovariance.push_back(0.0);
//                F >> inputFloatAmbCovariance[i];
//            }
//        } else if (n == 0) {
//            std::cerr << "No size given in the first line of input file.\n";
//        }
//    }
};


#endif //ILS_PPP_LAMBDATEST_H
