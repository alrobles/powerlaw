cmake_minimum_required(VERSION 3.21)
project(PowerLawFitterCpp VERSION 1.0 DESCRIPTION "A fitter for the discrete power-law distribution")
find_package(GSL REQUIRED)
find_package (Threads)

set(CMAKE_CXX_STANDARD 20)

option(CLI_BUILD "Build the CLI interface" OFF)
option(MATHLINK_BUILD "Build for Mathlink" OFF)
option(LIB_BUILD "Build a static library" OFF)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE EQUAL Release)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

add_executable(PowerLawFitterCppApp
            src/Zeta.h
            src/Zeta.cpp
            src/RandomGen.cpp
            src/ThreadPool.h
            src/TestStatistics.cpp
            src/DiscreteDistributions.cpp
            src/VectorUtilities.h
            src/ProgressBar.h
            src/ProgressBar.cpp
            include/DiscreteDistributions.h
            include/RandomGen.h
            include/TestStatistics.h
            cli/CLIMain.cpp
            cli/CsvParser.h
            cli/OptionParser.h)
target_link_libraries(PowerLawFitterCppApp GSL::gsl GSL::gslcblas ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS PowerLawFitterCppApp
    COMPONENT linapp
    RUNTIME DESTINATION "/home/"
    LIBRARY DESTINATION "/home/"
    DESTINATION "/home/"
)

SET(CPACK_GENERATOR "DEB")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Angel Robles") #required
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libgsl-dev")
INCLUDE(CPack)
