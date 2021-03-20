#include <iostream>
#include <sstream>
#include <fstream>
#include <istream>
#include "Matrix.h"
#include "Lambda.h"
#include "LambdaTest.h"

int LambdaTest::test() {
    const std::vector<double> inputAmbiguity{-9.75792, 22.1086, -1.98908, 3.36186, 23.2148, 7.75073};
    const std::vector<double> inputFloatAmbCovariance{
            0.0977961, 0.0161137, 0.0468261, 0.0320695, 0.080857, 0.0376408,
            0.0161137, 0.0208976, 0.0185378, 0.00290225, 0.0111409, 0.0247762,
            0.0468261, 0.0185378, 0.0435412, 0.0227732, 0.0383208, 0.0382978,
            0.0320695, 0.00290225, 0.0227732, 0.0161712, 0.0273471, 0.0154774,
            0.080857, 0.0111409, 0.0383208, 0.0273471, 0.0672121, 0.0294637,
            0.0376408, 0.0247762, 0.0382978, 0.0154774, 0.0294637, 0.0392536};
    Matrix<double> floatAmbiguity;
    floatAmbiguity.pushBackColumn(inputAmbiguity);

    Matrix<double> floatAmbCovariance =
            Matrix<double>::pushMatrixFromVector(inputFloatAmbCovariance, inputAmbiguity.size());
    std::cout << floatAmbCovariance;

    Lambda lambdaSolution;
    Matrix<int> integerAmbiguity(floatAmbiguity.getRows(), 1);
    integerAmbiguity =
            lambdaSolution.computeIntegerSolution(floatAmbiguity, floatAmbCovariance);

    std::cout << "Float Ambiguity: \n" << floatAmbiguity.transpose();
    std::cout << "Integer Ambiguity: \n" << integerAmbiguity.transpose();

    return 0;
}

int LambdaTest::testFromFile() {
    std::ifstream inputFile;
    std::vector<double> inputAmbiguity;
    std::vector<double> inputFloatAmbCovariance;
    int n{0};
    inputFile.open("input.txt");
    inputFile >> n;
    if (n != 0) {
        for (int i = 0; i < n; ++i) {
            inputAmbiguity.push_back(0.0);
            inputFile >> inputAmbiguity[i];
        }
        for (int i = 0; i < n * n; ++i) {
            inputFloatAmbCovariance.push_back(0.0);
            inputFile >> inputFloatAmbCovariance[i];
        }
    } else if (n == 0) {
        std::cerr << "No size given in the first line of input file.\n";
    }
    inputFile.close();
    Matrix<double> floatAmbiguity;
    floatAmbiguity.pushBackColumn(inputAmbiguity);
    std::cout << floatAmbiguity << '\n';

    Matrix<double> floatAmbCovariance =
            Matrix<double>::pushMatrixFromVector(inputFloatAmbCovariance, inputAmbiguity.size());
    std::cout << floatAmbCovariance;

    Lambda lambdaSolution;
    Matrix<int> integerAmbiguity(floatAmbiguity.getRows(), 1);
    integerAmbiguity =
            lambdaSolution.computeIntegerSolution(floatAmbiguity, floatAmbCovariance);

    std::cout << "Float Ambiguity: \n" << floatAmbiguity.transpose();
    std::cout << "Integer Ambiguity: \n" << integerAmbiguity.transpose();
    return 0;
}