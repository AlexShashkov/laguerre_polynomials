cmake_minimum_required(VERSION 3.26)
project(laguerre_polynomials)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(headers)

add_executable(test_roots_poly
        headers/ExtendedFMA.h
        headers/FastFourier.h
        headers/Polynomial.h
        roots.poly.h
        test_roots_poly.cpp)

add_executable(test_fpml
        headers/ExtendedFMA.h
        headers/FastFourier.h
        headers/Polynomial.h
        fpml.h
        test_fpml.cpp)

add_executable(laguerre_polynomials
        headers/ExtendedFMA.h
        headers/FastFourier.h
        headers/Polynomial.h
        roots.poly.h
        test_fpml.cpp)
