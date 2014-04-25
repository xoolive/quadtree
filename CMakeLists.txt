cmake_minimum_required(VERSION 2.6)
project (detection CXX)

option (BUILD_TESTS "Build tests" ON)

find_program (INTEL_COMPILER icc)
if (INTEL_COMPILER)
  include (CMakeForceCompiler)
  cmake_force_c_compiler (icc "Intel C Compiler")
  cmake_force_cxx_compiler (icpc "Intel C++ Compiler")
endif (INTEL_COMPILER)

add_library (quadtree STATIC
  neighbour.cpp
  quadtree.cpp)

enable_testing ()

if (BUILD_TESTS)
  add_subdirectory (tests)
endif (BUILD_TESTS)

