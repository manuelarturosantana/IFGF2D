cmake_minimum_required(VERSION 3.16)


set(COMMON_COMPILE_FLAGS "-O3" "-Wall" "-Wextra")


# set the project name and version
project(LayerPotentials VERSION 1.0)

message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")

#/usr/bin/c++
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
enable_language(Fortran) 
set(CMAKE_Fortran_COMPILER ~/intel/oneapi/compiler/2025.0/bin/ifx)
message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")

################################### SOURCE AND LINKER PATHS ###############################
set(COMP_BESSEL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/complex_bessel-master/build/lib64/libcomplex_bessel.so")

set(GEOMETRY_SOURCE_DIR Geometry)
set(GEOMETRY_SOURCES 
    ${GEOMETRY_SOURCE_DIR}/ClosedCurve.cpp)

set(FORWARD_MAP_SOURCE_DIR ForwardMap)
set(FORWARD_MAP_SOURCES
    ${FORWARD_MAP_SOURCE_DIR}/LinearOperator.cpp)

# $ENV{MKLROOT}/include
# $ENV{I_MPI_ROOT}/include
include_directories(/opt/Eigen/)
# set(COMMON_INCLUDE_DIRS
#     /opt/Eigen/
# )
##########################################################################################

# add_executable(IFGF_main
#                 IFGF_Source/main.cpp)


# target_link_libraries(IFGF_main PRIVATE ${COMP_BESSEL_PATH})
# target_compile_options(IFGF_main PRIVATE ${COMMON_COMPILE_FLAGS})

# add_executable(IgnoreSingNSing
#                 ${GEOMETRY_SOURCES} tests/IFGF/IgnoreSingNSing.cpp)
# target_link_libraries(IgnoreSingNSing PRIVATE ${COMP_BESSEL_PATH})
# target_compile_options(IgnoreSingNSing PRIVATE ${COMMON_COMPILE_FLAGS})




##########################################################################################
################################### TESTS ################################################
##########################################################################################


############################### FORWARD MAP TESTS #######################################
# add_executable(InitForwardMap 
#                 ${GEOMETRY_SOURCES} tests/ForwardMap/InitializeForwardMap.cpp)
# target_compile_options(InitForwardMap PRIVATE ${COMMON_COMPILE_FLAGS})

# add_executable(SinglePointSingCoeffs
#                 ${GEOMETRY_SOURCES} tests/ForwardMap/SinglePointSingCoeffs.cpp)
# target_compile_options(SinglePointSingCoeffs PRIVATE ${COMMON_COMPILE_FLAGS})
# target_link_libraries(SinglePointSingCoeffs PRIVATE ${COMP_BESSEL_PATH})

# add_executable(ComputeAx
#                 ${GEOMETRY_SOURCES} tests/ForwardMap/ComputeAx.cpp)
# target_compile_options(ComputeAx PRIVATE ${COMMON_COMPILE_FLAGS})
# target_link_libraries(ComputeAx PRIVATE ${COMP_BESSEL_PATH})

add_executable(SolveAndPropagate
                ${GEOMETRY_SOURCES} ${FORWARD_MAP_SOURCES} tests/ForwardMap/SolveAndPropagate.cpp)
target_compile_options(SolveAndPropagate PRIVATE ${COMMON_COMPILE_FLAGS})
target_link_libraries(SolveAndPropagate PRIVATE ${COMP_BESSEL_PATH})



################################# GEOMETRY TESTS #########################################
# add_executable(ClosedCurvePatchSplitting
#                 ${GEOMETRY_SOURCES} tests/Geometry/ClosedCurvePatchSplitting.cpp)
# target_compile_options(ClosedCurvePatchSplitting PRIVATE ${COMMON_COMPILE_FLAGS})

# add_executable(BoundingBoxTests
#                 ${GEOMETRY_SOURCES} tests/Geometry/BoundingBoxTests.cpp)
# target_compile_options(BoundingBoxTests PRIVATE ${COMMON_COMPILE_FLAGS})

# add_executable(ComputeBoundingBoxTests
#                 ${GEOMETRY_SOURCES} tests/Geometry/ComputeBoundingBoxTests.cpp)
# # target_include_directories(ComputeBoundingBoxTests PRIVATE ${COMMON_INCLUDE_DIRS})
# target_compile_options(ComputeBoundingBoxTests PRIVATE ${COMMON_COMPILE_FLAGS})

# ################################ UTILS TESTS #############################################
# add_executable(Cheb1d_test tests/utils/Cheb1d_test)
# target_include_directories(Cheb1d_test PRIVATE ${COMMON_INCLUDE_DIRS})
# target_compile_options(Cheb1d_test PRIVATE ${COMMON_COMPILE_FLAGS})

# add_executable(SingularIntegral tests/utils/SingularIntegral)
# target_compile_options(SingularIntegral PRIVATE ${COMMON_COMPILE_FLAGS})

