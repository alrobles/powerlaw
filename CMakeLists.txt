cmake_minimum_required(VERSION 3.21)
project(PowerLawFitterCpp VERSION 1.0 DESCRIPTION "A fitter for the discrete power-law distribution")
find_package(GSL REQUIRED)
find_package (Threads)

set(CMAKE_CXX_STANDARD 20)

option(CLI_BUILD "Build the CLI interface" OFF)
option(MATHLINK_BUILD "Build for Mathlink" OFF)
option(LIB_BUILD "Build a static library" OFF)

if(NOT CLI_BUILD AND NOT MATHLINK_BUILD AND NOT LIB_BUILD)
    set(CLI_BUILD ON)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE EQUAL Release)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

if (CLI_BUILD)
    add_executable(PowerLawFitterCpp
            src/Zeta.h
            src/Zeta.cpp
            src/RandomGen.cpp
            src/ThreadPool.h
            src/TestStatistics.cpp
            src/DiscreteDistributions.cpp
            src/VectorUtilities.h
            include/DiscreteDistributions.h
            include/RandomGen.h
            include/TestStatistics.h
            cli/CLIMain.cpp
            cli/CsvParser.h
            cli/OptionParser.h)

elseif(LIB_BUILD)
    add_library(PowerLawFitterCpp SHARED
            src/Zeta.h
            src/Zeta.cpp
            src/RandomGen.cpp
            src/ThreadPool.h
            src/TestStatistics.cpp
            src/DiscreteDistributions.cpp
            src/VectorUtilities.h
            include/DiscreteDistributions.h
            include/RandomGen.h
            include/TestStatistics.h)
    set_target_properties(PowerLawFitterCpp PROPERTIES VERSION ${PROJECT_VERSION})

elseif(MATHLINK_BUILD)
    set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/Mathematica" ${CMAKE_MODULE_PATH})
    find_package(Mathematica COMPONENTS MathLink)
    if (Mathematica_MathLink_DEFINITIONS)
        add_definitions(${Mathematica_MathLink_DEFINITIONS})
    endif()
    include_directories(${Mathematica_INCLUDE_DIRS})

    Mathematica_MathLink_MPREP_TARGET(mathlink/Mathlink.tm OUTPUT mathlink/Mathlinktm.cpp)

    add_executable(PowerLawFitterCpp
            src/Zeta.h
            src/Zeta.cpp
            src/RandomGen.cpp
            src/ThreadPool.h
            src/TestStatistics.cpp
            src/DiscreteDistributions.cpp
            src/VectorUtilities.h
            include/DiscreteDistributions.h
            include/RandomGen.h
            include/TestStatistics.h
            mathlink/MathlinkInterface.cpp
            mathlink/Mathlinktm.cpp
            )

    target_link_libraries(PowerLawFitterCpp ${Mathematica_MathLink_LIBRARIES})
endif()

target_link_libraries(PowerLawFitterCpp GSL::gsl GSL::gslcblas ${CMAKE_THREAD_LIBS_INIT})