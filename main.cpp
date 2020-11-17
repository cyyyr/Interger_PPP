#include <iostream>
#include <fstream>
#include <sstream>
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
    Matrix<double> floatAmbiguity(6,1);
    Matrix<double> floatAmbCovariance(6,6);
    std::stringstream ss;
    ss << -9.75792 << 22.1086 << -1.98908 << 3.36186
       <<  23.2148 << 7.75073;
    ss >> floatAmbiguity;

    ss << 0.0977961 << 0.0161137  << 0.0468261 << 0.0320695  << 0.080857  << 0.0376408
       << 0.0161137 << 0.0208976  << 0.0185378 << 0.00290225 << 0.0111409 << 0.0247762
       << 0.0468261 << 0.0185378  << 0.0435412 << 0.0227732  << 0.0383208 << 0.0382978
       << 0.0320695 << 0.00290225 << 0.0227732 << 0.0161712  << 0.0273471 << 0.0154774
       << 0.080857  << 0.0111409  << 0.0383208 << 0.0273471  << 0.0672121 << 0.0294637
       << 0.0376408 << 0.0247762  << 0.0382978 << 0.0154774  << 0.0294637 << 0.0392536;

    Lambda LambdaResolution;
    Matrix<int> integerAmbiguity(floatAmbiguity.getRows(), 0);
    integerAmbiguity =
            LambdaResolution.resolveIntegerAmbiguity(floatAmbiguity, floatAmbCovariance);

    std::cout << "Float Ambiguity: \n" << floatAmbiguity.transpose() << '\n';
    std::cout << "Integer Ambiguity: \n" << integerAmbiguity.transpose() << '\n';

//    std::cout << mat(0,0) << ' '
//              << mat(1,0) << ' '
//              << mat(2,0) << ' '
//              << mat(3,0) << ' '
//              << mat(4,0) << ' '
//              << mat(5,0) << '\n';

    std::cin.get();
}
