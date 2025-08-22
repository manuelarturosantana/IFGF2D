cmake_minimum_required(VERSION 3.16)


set(COMMON_COMPILE_FLAGS "-O3")


# set the project name and version
project(LayerPotentials VERSION 1.0)

message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")

#/usr/bin/c++
# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
enable_language(Fortran) 
set(CMAKE_Fortran_COMPILER ~/intel/oneapi/compiler/2025.0/bin/ifx)
message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")

################################### SOURCE AND LINKER PATHS ###############################
set(COMP_BESS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/complex_bessel-master/build/lib64/libcomplex_bessel.so")

set(GEOMETRY_SOURCE_DIR Geometry)
set(GEOMETRY_SOURCES 
    ${GEOMETRY_SOURCE_DIR}/ClosedCurve.cpp
    ${GEOMETRY_SOURCE_DIR}/Patch.cpp)

# $ENV{MKLROOT}/include
# $ENV{I_MPI_ROOT}/include
include_directories(/opt/Eigen/)
# set(COMMON_INCLUDE_DIRS
#     /opt/Eigen/
# )
##########################################################################################

add_executable(IFGF_main
                IFGF_Source/main.cpp)


target_link_libraries(IFGF_main PRIVATE ${COMP_BESS_PATH})
target_compile_options(IFGF_main PRIVATE -Wall -Wextra)




##########################################################################################
################################### TESTS ################################################
##########################################################################################

################################# GEOMETRY TESTS #########################################
add_executable(ClosedCurvePatchSplitting
                ${GEOMETRY_SOURCES} tests/Geometry/ClosedCurvePatchSplitting.cpp)
target_compile_options(ClosedCurvePatchSplitting PRIVATE ${COMMON_COMPILE_FLAGS})

add_executable(BoundingBoxTests
                ${GEOMETRY_SOURCES} tests/Geometry/BoundingBoxTests.cpp)
target_compile_options(BoundingBoxTests PRIVATE ${COMMON_COMPILE_FLAGS})

add_executable(ComputeBoundingBoxTests
                ${GEOMETRY_SOURCES} tests/Geometry/ComputeBoundingBoxTests.cpp)
# target_include_directories(ComputeBoundingBoxTests PRIVATE ${COMMON_INCLUDE_DIRS})
target_compile_options(ComputeBoundingBoxTests PRIVATE ${COMMON_COMPILE_FLAGS})

################################ UTILS TESTS #############################################
add_executable(Cheb1d_test tests/utils/Cheb1d_test)
# target_include_directories(Cheb1d_test PRIVATE ${COMMON_INCLUDE_DIRS})
target_compile_options(Cheb1d_test PRIVATE ${COMMON_COMPILE_FLAGS})