#include <iostream>
#include <sstream>
#include "Matrix.h"
#include "Lambda.h"
#include "LambdaTest.h"

int LambdaTest::test(){
    try {
        Matrix<double> floatAmbiguity;
        floatAmbiguity.pushBackColumn({-9.75792, 22.1086, -1.98908, 3.36186, 23.2148, 7.75073});

        Matrix<double> floatAmbCovariance;
        floatAmbCovariance.pushBackRow({0.0977961, 0.0161137, 0.0468261, 0.0320695, 0.080857, 0.0376408});
        floatAmbCovariance.pushBackRow({0.0161137, 0.0208976, 0.0185378, 0.00290225, 0.0111409, 0.0247762});
        floatAmbCovariance.pushBackRow({0.0468261, 0.0185378, 0.0435412, 0.0227732, 0.0383208, 0.0382978});
        floatAmbCovariance.pushBackRow({0.0320695, 0.00290225, 0.0227732, 0.0161712, 0.0273471, 0.0154774});
        floatAmbCovariance.pushBackRow({0.080857, 0.0111409, 0.0383208, 0.0273471, 0.0672121, 0.0294637});
        floatAmbCovariance.pushBackRow({0.0376408, 0.0247762, 0.0382978, 0.0154774, 0.0294637, 0.0392536});

        Lambda lambdaSolution;
        Matrix<int> integerAmbiguity(floatAmbiguity.getRows(), 1);
        integerAmbiguity =
                lambdaSolution.computeIntegerSolution(floatAmbiguity, floatAmbCovariance);

        std::cout << "Float Ambiguity: \n" << floatAmbiguity.transpose();
        std::cout << "Integer Ambiguity: \n" << integerAmbiguity.transpose();

        return 0;
    } catch (...) {
        std::cerr << "Testing error!\n";
        return -1;
    }
}