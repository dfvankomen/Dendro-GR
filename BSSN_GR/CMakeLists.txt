#
# CMAKE CONFIGURATION FILE FOR BSSN SOLVER
#
# This file in particular contains all of the options that apply specificially
# to the BSSN Solver binary as it is compiled.
#
# Note that all options and compiler flags set in the ../CMakeLists.txt file
# (the one in the root directory of the project) *are carried through* here.
# Meaning, at this point, the C++ standard is set, OpenMP, MPI, GSL, and other
# libraries are found and linked. Also, general Dendro options, such as ordering
# options, dimensionality, CMake's build type, and others will also propagate
# through to here.
#
#

# The BSSN Solver requires TOML in order to read from parameters (as an option)
# the toml11 library handles the reading for us.
# set(EXTERNAL_LIB_HEADERS ${CMAKE_SOURCE_DIR}/lib/toml11)

# these were already found, but it's good to repeat them
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(GSL REQUIRED)

# only BSSN (GR) options
option(BSSN_PROFILE_HUMAN_READABLE "enables human readable profile outputs" ON)
option(BSSN_CONSEC_COMM_SELECT "rank selection performed using consecative ranks (works with arbitary number of procs.) if off ranks selected in binary tree fashion. (works only with powers of 2)" ON)
option(BSSN_ENABLE_VTU_OUTPUT "Enables vtu output" ON)
option(BSSN_COMPUTE_CONSTRAINTS "Enable constraint equation computations" ON)
option(BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT "write constrints to vtu" ON)
option(BSSN_ETA_FUNCTION "Use function as ETA damping" OFF)
option(BSSN_KERR_SCHILD_TEST "compare with the analytical Kerr-Schild sol. for a single BH" OFF)
option(BSSN_ENABLE_CUDA "enable rhs computation with GPUs" OFF)
#option(BSSN_RHS_STAGED_COMP "enable rhs staged computation" OFF)
option(BSSN_GAUGE_ROCHESTER "enable rochester gauge for bssn equations" OFF)
option(BSSN_EXTRACT_BH_LOCATIONS "compute the BH locations as time evolves" ON)
option(BSSN_REFINE_BASE_EH "enable refine based on BH event horizons" OFF)
option(BSSN_EXTRACT_GRAVITATIONAL_WAVES "extract GW waves" ON)
option(BSSN_VERBOSE_LOGGING "Enable more verbose logging (only useful during debug)" OFF)
option(BSSN_ENABLE_SSL_HD "Enable SSL and HD formalism modifications" ON)

option(BSSN_USE_4TH_ORDER_DERIVS "use 4th order derivs" OFF)
option(BSSN_USE_6TH_ORDER_DERIVS "use 6th order derivs" ON)
option(BSSN_USE_8TH_ORDER_DERIVS "use 8th order derivs" OFF)

option(BSSN_ENABLE_AVX "use vectorized computation for RHS (intel compilers only)" OFF)


if(BSSN_ENABLE_AVX)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_definitions(-DBSSN_ENABLE_AVX)
    add_definitions(-DDERIV_ENABLE_AVX)
  else()
    message("AVX vectorization currently supported to intel compilers only. use export CC=icc and export CXX=icpc to use intel compilers. AVX vectorization not enabled. ")
  endif()

endif(BSSN_ENABLE_AVX)

if(BSSN_USE_4TH_ORDER_DERIVS)
  set(BSSN_USE_6TH_ORDER_DERIVS OFF)
  set(BSSN_USE_8TH_ORDER_DERIVS OFF)
elseif(BSSN_USE_6TH_ORDER_DERIVS)
  set(BSSN_USE_4TH_ORDER_DERIVS OFF)
  set(BSSN_USE_8TH_ORDER_DERIVS OFF)
elseif(BSSN_USE_8TH_ORDER_DERIVS)
  set(BSSN_USE_4TH_ORDER_DERIVS OFF)
  set(BSSN_USE_6TH_ORDER_DERIVS OFF)
endif()

set(BSSN_SRC_FILES
  src/parameters.cpp
  src/grUtils.cpp
  src/rhs.cpp
  src/derivs.cpp
  src/physcon.cpp
  src/profile_params.cpp
  src/bssn_constraints.cpp
  src/TPCoordTransf.cpp
  src/TPEquations.cpp
  src/TPFuncAndJacobian.cpp
  src/TPNewton.cpp
  src/TPUtilities.cpp
  src/TwoPunctures.cpp
  src/dataUtils.cpp
  src/hadrhs.F90
  src/trumpet.cpp
)

add_library(bssn_common OBJECT ${BSSN_SRC_FILES})

set(BSSN_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts
)

if(DEFINED CHECK_GIT_INCLUDE_PATH)
  list(APPEND BSSN_INCLUDE_DIRS ${CHECK_GIT_INCLUDE_PATH})
endif()

target_include_directories(bssn_common PUBLIC ${BSSN_INCLUDE_DIRS})
target_compile_options(bssn_common PUBLIC -fext-numeric-literals)
target_link_libraries(bssn_common PUBLIC dendro5 MPI::MPI_CXX GSL::gsl toml11::toml11 dendro_git_version_and_date)

macro(bssn_add_def OPTION_NAME DEF_NAME)
  if(${OPTION_NAME})
    target_compile_definitions(bssn_common PUBLIC ${DEF_NAME})
  endif()
endmacro()

bssn_add_def(BSSN_USE_4TH_ORDER_DERIVS BSSN_USE_4TH_ORDER_DERIVS)
bssn_add_def(BSSN_USE_6TH_ORDER_DERIVS BSSN_USE_6TH_ORDER_DERIVS)
bssn_add_def(BSSN_USE_8TH_ORDER_DERIVS BSSN_USE_8TH_ORDER_DERIVS)
bssn_add_def(BSSN_PROFILE_HUMAN_READABLE BSSN_PROFILE_HUMAN_READABLE)
bssn_add_def(BSSN_CONSEC_COMM_SELECT BSSN_CONSEC_COMM_SELECT)
bssn_add_def(BSSN_ENABLE_VTU_OUTPUT BSSN_ENABLE_VTU_OUTPUT)
bssn_add_def(BSSN_COMPUTE_CONSTRAINTS BSSN_COMPUTE_CONSTRAINTS)
bssn_add_def(BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT)
bssn_add_def(BSSN_ETA_FUNCTION USE_ETA_FUNC)
bssn_add_def(BSSN_KERR_SCHILD_TEST BSSN_KERR_SCHILD_TEST)
bssn_add_def(BSSN_ENABLE_CUDA BSSN_ENABLE_CUDA)
# bssn_add_def(BSSN_RHS_STAGED_COMP BSSN_RHS_STAGED_COMP)
bssn_add_def(BSSN_EXTRACT_BH_LOCATIONS BSSN_EXTRACT_BH_LOCATIONS)
bssn_add_def(BSSN_REFINE_BASE_EH BSSN_REFINE_BASE_EH)
bssn_add_def(BSSN_GAUGE_ROCHESTER USE_ROCHESTER_GAUGE)
bssn_add_def(BSSN_VERBOSE_LOGGING BSSN_VERBOSE_LOGGING)
bssn_add_def(BSSN_ENABLE_SSL_HD BSSN_ENABLE_SSL_HD)
bssn_add_def(BSSN_EXTRACT_GRAVITATIONAL_WAVES BSSN_EXTRACT_GRAVITATIONAL_WAVES)

if(BSSN_ENABLE_AVX AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  target_compile_definitions(bssn_common PUBLIC BSSN_ENABLE_AVX DERIV_ENABLE_AVX)
endif()


# CUDA setup:
if(WITH_CUDA)
  set(CUDA_GR_INC
      cuda/include/test_utils.h
      cuda/include/bssn_kernels.cuh
      cuda/include/bssnrhs_evar_derivs.cuh)

  cuda_add_library(bssn_cuda ${BSSN_SRC_FILES} ${CUDA_GR_INC})
  target_include_directories(bssn_cuda PRIVATE
    ${BSSN_INCLUDE_DIRS}
    ${DENDRO_dendrolib_DIR}/include
    ${DENDRO_dendrolib_DIR}/GPU/include
    ${MPI_INCLUDE_PATH}
    ${GSL_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include
  )
  target_link_libraries(bssn_cuda toml11::toml11)
endif()

### EXECUTABLES:
add_executable(bssnSolver
    src/bssngr_main.cpp
    src/bssnAEH.cpp
    src/bssnCtx.cpp
)
target_link_libraries(bssnSolver PRIVATE bssn_common)
target_link_libraries(bssnSolver PRIVATE dendro5 MPI::MPI_CXX GSL::gsl toml11::toml11 m dendro_git_version_and_date)
target_compile_options(bssnSolver PRIVATE -fext-numeric-literals)

if(WITH_CUDA)
  target_include_directories(bssnSolver PRIVATE ${CUDA_INCLUDE_DIRS})
  target_link_libraries(bssnSolver PRIVATE bssn_cuda ${CUDA_LIBRARIES})
endif()

# --- TPID ---
add_executable(tpid src/TPID.cpp)
target_link_libraries(tpid PRIVATE bssn_common)
target_link_libraries(tpid PRIVATE dendro5 MPI::MPI_CXX GSL::gsl toml11::toml11 m dendro_git_version_and_date)
target_compile_options(tpid PRIVATE -fext-numeric-literals)
if(WITH_CUDA)
  target_include_directories(tpid PRIVATE ${CUDA_INCLUDE_DIRS})
  target_link_libraries(tpid PRIVATE bssn_cuda ${CUDA_LIBRARIES})
endif()

# TESTS AND OTHER EXECUTABLES
add_executable(bssnWAMRTest src/gr_wamr_test.cpp)
target_link_libraries(bssnWAMRTest PRIVATE bssn_common)
target_link_libraries(bssnWAMRTest PRIVATE dendro5 MPI::MPI_CXX GSL::gsl toml11::toml11 m dendro_git_version_and_date)
target_compile_options(bssnWAMRTest PRIVATE -fext-numeric-literals)
if(WITH_CUDA)
  target_include_directories(bssnWAMRTest PRIVATE ${CUDA_INCLUDE_DIRS})
  target_link_libraries(bssnWAMRTest PRIVATE bssn_cuda ${CUDA_LIBRARIES})
endif()


add_executable(checkDerivs src/checkDerivatives.cpp)
target_link_libraries(checkDerivs PRIVATE bssn_common)
target_link_libraries(checkDerivs PRIVATE dendro5 MPI::MPI_CXX GSL::gsl toml11::toml11 m dendro_git_version_and_date)
target_compile_options(checkDerivs PRIVATE -fext-numeric-literals)
if(WITH_CUDA)
  target_include_directories(checkDerivs PRIVATE ${CUDA_INCLUDE_DIRS})
  target_link_libraries(checkDerivs PRIVATE bssn_cuda ${CUDA_LIBRARIES})
endif()


# CUDA ONLY
if(WITH_CUDA)
  cuda_add_executable(derivs_test cuda/src/derivs_test.cu)
  target_include_directories(derivs_test PRIVATE ${BSSN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
  target_link_libraries(derivs_test bssn_cuda dendro5_cuda ${CUDA_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m toml11::toml11)

  cuda_add_executable(rhsgpu_test cuda/src/rhs_test.cu)
  target_include_directories(rhsgpu_test PRIVATE ${BSSN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
  target_link_libraries(rhsgpu_test bssn_cuda dendro5_cuda ${CUDA_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m toml11::toml11)

  add_executable(rhscpu_test cuda/src/rhs_test.cpp)
  target_link_libraries(rhscpu_test PRIVATE bssn_common)
  target_link_libraries(rhscpu_test PRIVATE dendro5 MPI::MPI_CXX GSL::gsl toml11::toml11 m dendro_git_version_and_date)
  target_compile_options(rhscpu_test PRIVATE -fext-numeric-literals)

  cuda_add_executable(bssnSolverCUDA cuda/src/gr_cuda.cu cuda/src/bssnCtxGPU.cu)
  target_include_directories(bssnSolverCUDA PRIVATE ${BSSN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
  target_link_libraries(bssnSolverCUDA bssn_cuda dendro5_cuda ${CUDA_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m toml11::toml11)

  cuda_add_executable(bssnWSTestCUDA cuda/src/gr_scaling.cu cuda/src/bssnCtxGPU.cu)
  target_include_directories(bssnWSTestCUDA PRIVATE ${BSSN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
  target_link_libraries(bssnWSTestCUDA bssn_cuda dendro5_cuda ${CUDA_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m toml11::toml11)
endif()







# set(SOURCE_FILES src/gr_scaling.cpp include/gr.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(bssnScalingTest ${SOURCE_FILES})
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/include)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/test/include)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/examples/include)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/FEM/include)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/ODE/include)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/LinAlg/include)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/IO/vtk/include)
# target_include_directories(bssnScalingTest PRIVATE ${DENDRO_dendrolib_DIR}/IO/zlib/inc)
# target_include_directories(bssnScalingTest PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(bssnScalingTest PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#     target_include_directories(bssnScalingTest PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(bssnScalingTest bssn_cuda)
#     target_link_libraries(bssnScalingTest ${CUDA_LIBRARIES})
# endif()
# target_link_libraries(bssnScalingTest dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


# set(SOURCE_FILES src/gr_nuts.cpp include/gr.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(bssnSolverNUTS ${SOURCE_FILES})

# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/test/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/examples/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/FEM/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/ODE/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/LinAlg/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/IO/vtk/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${DENDRO_dendrolib_DIR}/IO/zlib/inc)
# target_include_directories(bssnSolverNUTS PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(bssnSolverNUTS PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#     target_include_directories(bssnSolverNUTS PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(bssnSolverNUTS bssn_cuda)
#     target_link_libraries(bssnSolverNUTS ${CUDA_LIBRARIES})
# endif()

# target_link_libraries(bssnSolverNUTS dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


# if(WITH_CUDA)
#     set(SOURCE_FILES cuda_gr/src/gpuBSSNExample.cpp ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#     add_executable(cudaBSSN ${SOURCE_FILES})
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/include/cuda)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/test/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/examples/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/FEM/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/ODE/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/LinAlg/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/IO/vtk/include)
#     target_include_directories(cudaBSSN PRIVATE ${DENDRO_dendrolib_DIR}/IO/zlib/inc)
#     target_include_directories(cudaBSSN PRIVATE ${MPI_INCLUDE_PATH})
#     target_include_directories(cudaBSSN PRIVATE ${GSL_INCLUDE_DIRS})
#     target_include_directories(cudaBSSN PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(cudaBSSN bssn_cuda)
#     target_link_libraries(cudaBSSN ${CUDA_LIBRARIES})
#     target_link_libraries(cudaBSSN dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES}  m)
#     add_dependencies(cudaBSSN BSSN_CODEGEN)
#     #add_dependencies(cudaBSSN GW_EXTRACT)
# endif()



# if(WITH_CUDA)
# set(SOURCE_FILES cuda_gr/src/gpuTestBed.cpp cuda_gr/include/gpuTestBed.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(gpuTestBed ${SOURCE_FILES})
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/include/cuda)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/test/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/examples/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/FEM/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/ODE/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/LinAlg/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/IO/vtk/include)
#     target_include_directories(gpuTestBed PRIVATE ${DENDRO_dendrolib_DIR}/IO/zlib/inc)
#     target_include_directories(gpuTestBed PRIVATE ${MPI_INCLUDE_PATH})
#     target_include_directories(gpuTestBed PRIVATE ${GSL_INCLUDE_DIRS})
#     target_include_directories(gpuTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(gpuTestBed bssn_cuda)
#     target_link_libraries(gpuTestBed ${CUDA_LIBRARIES})
#     target_link_libraries(gpuTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)
#     add_dependencies(gpuTestBed BSSN_CODEGEN)
#     #add_dependencies(gpuTestBed GW_EXTRACT)
# endif()


# set(SOURCE_FILES src/rhsTestBed.cpp include/rhsTestBed.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(rhsTestBed ${SOURCE_FILES})
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/include/cuda)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/test/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/examples/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/FEM/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/ODE/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/LinAlg/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/IO/vtk/include)
# target_include_directories(rhsTestBed PRIVATE ${DENDRO_dendrolib_DIR}/IO/zlib/inc)
# target_include_directories(rhsTestBed PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(rhsTestBed PRIVATE ${GSL_INCLUDE_DIRS})

# if(WITH_CUDA)
#     target_include_directories(rhsTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(rhsTestBed bssn_cuda)
#     target_link_libraries(rhsTestBed ${CUDA_LIBRARIES})
# endif()

# target_link_libraries(rhsTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


# NOTE: APRIL 10 DFVK - I turned this off for now due to some BSSN tests with extra values

# set(SOURCE_FILES src/zipUnzipTest.cpp include/zipUnzipTest.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(zipUnzipTest ${SOURCE_FILES})
# target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/include)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/test/include)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/examples/include)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/FEM/include)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/ODE/include)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/LinAlg/include)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/IO/vtk/include)
# target_include_directories(zipUnzipTest PRIVATE ${DENDRO_dendrolib_DIR}/IO/zlib/inc)
# target_include_directories(zipUnzipTest PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(zipUnzipTest PRIVATE ${GSL_INCLUDE_DIRS})
# target_include_directories(zipUnzipTest PRIVATE ${EXTERNAL_LIB_HEADERS})  # includes for the TOML11 library (and others)
#
# if(WITH_CUDA)
#     target_include_directories(zipUnzipTest PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(zipUnzipTest bssn_cuda)
#     target_link_libraries(zipUnzipTest ${CUDA_LIBRARIES})
# endif()
#
# target_link_libraries(zipUnzipTest dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m dendro_git_version_and_date)






