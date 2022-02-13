cmake_minimum_required(VERSION 3.21)
project(PowerLawFitterCpp)
find_package(GSL REQUIRED)
find_package (Threads)

set(CMAKE_CXX_STANDARD 20)

add_executable(PowerLawFitterCpp
        Main.cpp
        Zeta.h
        Zeta.cpp
        RandomGen.cpp
        RandomGen.h
        ThreadPool.h
        TestStatistics.cpp
        TestStatistics.h
        DiscreteDistributions.cpp
        DiscreteDistributions.h)

target_link_libraries(PowerLawFitterCpp GSL::gsl GSL::gslcblas ${CMAKE_THREAD_LIBS_INIT})