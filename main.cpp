#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Matrix.h"
#include "Lambda.h"

int main() {
    //TestILSMethod();
    /* 1) Чтение из файла (см. Coursera)
     * 2) МНК
     * 3) Исправление данных
     * 4) Тест система
     * 5) График
    */

    Matrix<double> floatAmbiguity;
    floatAmbiguity.pushBackColumn({-9.75792, 22.1086, -1.98908, 3.36186, 23.2148, 7.75073});

    Matrix<double> floatAmbCovariance;
    floatAmbCovariance.pushBackRow({0.0977961, 0.0161137, 0.0468261, 0.0320695, 0.080857, 0.0376408});
    floatAmbCovariance.pushBackRow({0.0161137, 0.0208976, 0.0185378, 0.00290225, 0.0111409, 0.0247762});
    floatAmbCovariance.pushBackRow({0.0468261, 0.0185378, 0.0435412, 0.0227732, 0.0383208, 0.0382978});
    floatAmbCovariance.pushBackRow({0.0320695, 0.00290225, 0.0227732, 0.0161712, 0.0273471, 0.0154774});
    floatAmbCovariance.pushBackRow({0.080857, 0.0111409, 0.0383208, 0.0273471, 0.0672121, 0.0294637});
    floatAmbCovariance.pushBackRow({0.0376408, 0.0247762, 0.0382978, 0.0154774, 0.0294637, 0.0392536});


    Lambda LambdaResolution;
    Matrix<int> integerAmbiguity(floatAmbiguity.getRows(), 1);
    integerAmbiguity =
            LambdaResolution.resolveIntegerAmbiguity(floatAmbiguity, floatAmbCovariance);

    std::cout << "Float Ambiguity: \n" << floatAmbiguity.transpose() << '\n';
    std::cout << "Integer Ambiguity: \n" << integerAmbiguity.transpose() << '\n';

//    std::cout << floatAmbiguity(0,0) << ' ' << floatAmbiguity(0,1) << '\n'
//              << floatAmbiguity(1,0) << ' ' << floatAmbiguity(1,1) << '\n'
//              << floatAmbiguity(2,0) << ' ' << floatAmbiguity(2,1) << '\n'
//              << floatAmbiguity(3,0) << ' ' << floatAmbiguity(3,1) << '\n'
//              << floatAmbiguity(4,0) << ' ' << floatAmbiguity(4,1) << '\n'
//              << floatAmbiguity(5,0) << ' ' << floatAmbiguity(5,1) << '\n';
//
//    std::cout << "\n\n\n";
//
//    std::cout << floatAmbCovariance(0,0) << ' ' << floatAmbCovariance(0,1) << ' ' << floatAmbCovariance(0,2) << ' '
//              << floatAmbCovariance(0,3) << ' ' << floatAmbCovariance(0,4) << ' ' << floatAmbCovariance(0,5) << '\n'
//            << floatAmbCovariance(5,0) << ' ' << floatAmbCovariance(5,1) << ' ' << floatAmbCovariance(5,2) << ' '
//            << floatAmbCovariance(5,3) << ' ' << floatAmbCovariance(5,4) << ' ' << floatAmbCovariance(5,5) << '\n';

    //std::cin.get();
}
