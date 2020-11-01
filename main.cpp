#include <iostream>
#include <fstream>
#include "Matrix.h"

int main() {
    //TestILSMethod();
    /* 1) Чтение из файла (см. Coursera)
     * 2) МНК
     * 3) Исправление данных
     * 4) Тест система
     * 5) График
    */
    int n = 2, m = 2;
    Matrix<int> mat(2,2);
    std::cout << mat.getRows() << ',' << mat.getCols() << std::endl;

    std::cin.get();
}
