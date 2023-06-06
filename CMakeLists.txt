cmake_minimum_required(VERSION 3.12)

project("SSH-Hubbard"
  DESCRIPTION "SSH-Hubbard model code."
  LANGUAGES CXX
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")


# Build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(${CMAKE_BUILD_TYPE} MATCHES Debug)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

option(GQTEN_TIMING_MODE "Timing mode." OFF)
if(GQTEN_TIMING_MODE)
  add_definitions(-DGQTEN_TIMING_MODE)
endif()

option(GQMPS2_TIMING_MODE "Timing mode." OFF)
if(GQMPS2_TIMING_MODE)
  add_definitions(-DGQMPS2_TIMING_MODE)
endif()

option(GQTEN_MPI_TIMING_MODE "Timing mode." OFF)
if(GQTEN_MPI_TIMING_MODE)
  add_definitions(-DGQTEN_MPI_TIMING_MODE)
endif()

option(GQMPS2_MPI_TIMING_MODE "Timing mode." OFF)
if(GQMPS2_MPI_TIMING_MODE)
  add_definitions(-DGQMPS2_MPI_TIMING_MODE)
endif()

# Compilation and linking control.
set(CMAKE_CXX_STANDARD 17)
add_definitions(-Wall -g)
if(${CMAKE_BUILD_TYPE} MATCHES Release)
  add_compile_options(-O3)
endif()



if(APPLE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
    #Need test
    set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
  endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
    # May not work
    set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -liomp5 -lpthread -lm -ldl)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
    # Note as of Intel oneAPI 2021.2, "source /opt/intel/oneapi/compiler/latest/env/vars.sh"
    set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -Wl, -rpath $ENV{CMPLR_ROOT}/mac/compiler/lib/libiomp5.dylib -liomp5 -lpthread -lm -ldl)
  endif()
elseif(UNIX)
  set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "~/.local/")
  message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")

  #UNIX include APPLE, but we except it here
  if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
    set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
    # Link the Intel's OpenMP library to avoid performance issue when the library calls the MKL's gesdd function.
    set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -lpthread -lirc -lsvml -Wl,--no-as-needed -ldl)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
    set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lpthread -lm -ldl)
  endif()
endif()



find_path(TENSOR_HEADER_PATH "gqten")
find_path(MPS2_HEADER_PATH "gqmps2")
find_package(hptt)
find_package(Boost
  1.76.0
  REQUIRED 
  COMPONENTS serialization mpi)
find_package(MPI REQUIRED)
#set omp flag
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(OMP_FLAGS -qopenmp)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(OMP_FLAGS -fopenmp)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(OMP_FLAGS -fopenmp)
endif()

#find openmp?
add_compile_options( ${MATH_LIB_COMPILE_FLAGS} ${OMP_FLAGS})

include_directories(  ${MPS2_HEADER_PATH}
 ${TENSOR_HEADER_PATH}
 ${hptt_INCLUDE_DIR}
 ${Boost_INCLUDE_DIR}
 ${MPI_CXX_HEADER_DIR}
)

link_libraries(
  ${hptt_LIBRARY}
  ${MATH_LIB_LINK_FLAGS}
  Boost::serialization
  Boost::mpi
  ${MPI_CXX_LINK_FLAGS}
  ${MPI_mpi_LIBRARY}
)


add_executable(mpogen src/ssh2dmpogen.cpp src/operators.cpp)
add_executable(OBCmpogen src/OBCmpogen.cpp src/operators.cpp)
add_executable(holsteinmpogen src/holstein2dmpogen.cpp src/operators.cpp)
#add_executable(pla_mpogen pla_ssh2dmpogen.cpp operators.cpp)
#add_executable(oy_mpogen oy_ssh2dmpogen.cpp operators.cpp)
#add_executable(ox_mpogen ox_ssh2dmpogen.cpp operators.cpp)
#add_executable(pla_sshMC pla_sshMC.cpp operators.cpp myutil.cpp)
#add_executable(oy_sshMC oy_sshMC.cpp operators.cpp myutil.cpp)
#add_executable(ox_sshMC ox_sshMC.cpp operators.cpp myutil.cpp)
add_executable(OBCsshMC src/OBCsshMC.cpp src/operators.cpp src/myutil.cpp)
add_executable(vmps_ssh_pbc src/vmps_ssh_pbc.cpp src/operators.cpp src/myutil.cpp)
add_executable(holstein_parallel src/holstein_parallel.cpp src/operators.cpp src/myutil.cpp)
#add_executable(pla_parallel pla_parallel.cpp operators.cpp myutil.cpp)
#add_executable(oy_parallel oy_parallel.cpp operators.cpp myutil.cpp)
#add_executable(ox_parallel ox_parallel.cpp operators.cpp myutil.cpp)
add_executable(OBCparallel src/OBCparallel.cpp src/operators.cpp src/myutil.cpp)
add_executable(measure2 src/measure2.cpp src/operators.cpp src/myutil.cpp)
add_executable(OBCmeasure2 src/OBCmeasure2.cpp src/operators.cpp src/myutil.cpp)
add_executable(measure1 src/measure1.cpp src/operators.cpp src/myutil.cpp)
add_executable(OBCmeasure1 src/OBCmeasure1.cpp src/operators.cpp src/myutil.cpp)
add_executable(holstein_measure1 src/holstein_measure1.cpp src/operators.cpp src/myutil.cpp)
add_executable(holstein_measure2 src/holstein_measure2.cpp src/operators.cpp src/myutil.cpp)
add_executable(measureSC src/measureSC.cpp src/operators.cpp src/myutil.cpp)
add_executable(OBCmeasureSC src/OBCmeasureSC.cpp src/operators.cpp src/myutil.cpp)
add_executable(holstein_measureSC src/holstein_measureSC.cpp src/operators.cpp src/myutil.cpp)
#add_executable(OBCmeasureSC_ps OBCmeasureSC_plaquette_sym.cpp operators.cpp myutil.cpp)
#add_executable(oy_measureSC_ps oy_measureSC_pla.cpp operators.cpp myutil.cpp)
#add_executable(ox_measureSC_ps ox_measureSC_pla.cpp operators.cpp myutil.cpp)
#add_executable(ox_measureSC ox_measureSC.cpp operators.cpp myutil.cpp)
#add_executable(pla_measureSC_ps pla_measureSC_pla.cpp operators.cpp myutil.cpp)
add_executable(ssh2ds src/ssh2ds.cpp src/operators.cpp src/myutil.cpp)
add_executable(move_center src/move_center.cpp src/myutil.cpp)
add_executable(fix_mps src/fix_mps.cpp src/myutil.cpp)
add_executable(fix_mps2 src/fix_mps2.cpp src/myutil.cpp)
add_executable(fix_update_renv src/fix_update_renv.cpp src/myutil.cpp)
add_executable(fix_update_lenv src/fix_update_lenv.cpp src/myutil.cpp)
add_executable(continue src/continue.cpp src/operators.cpp src/myutil.cpp)
add_executable(mpi_test src/mpi_test.cpp)