cmake_minimum_required(VERSION 3.17)
project(ILS_PPP)

set(CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")


add_executable(ILS_PPP main.cpp Lambda.cpp Lambda.h Matrix.cpp Matrix.h LambdaTest.cpp LambdaTest.h)