# Compilation of a smart version of quadtrees specialised for tracking moving
# objects. When iterating over elements in the quadtree, a simple flag
# (boolean) indicates whether the element might have moved to a neighbouring
# subdivision.
#
# Xavier Olive, 28 nov. 2014

cmake_minimum_required (VERSION 2.8.12)
project (quadtree)

option (BUILD_TESTS "Build tests" ON)

find_program (INTEL_COMPILER icc)
if (INTEL_COMPILER)
  include (CMakeForceCompiler)
  cmake_force_c_compiler (icc "Intel C Compiler")
  cmake_force_cxx_compiler (icpc "Intel C++ Compiler")
endif (INTEL_COMPILER)

if (WIN32)
  if (MINGW)
    set (CMAKE_CXX_FLAGS "-std=c++0x")
  endif (MINGW)
else (WIN32)
  #set (CMAKE_CXX_FLAGS "-std=c++11")
  set (CMAKE_CXX_FLAGS "-std=c++0x -O2 -Wall -pedantic")
endif(WIN32)

add_library (smartquadtree STATIC
  neighbour.cpp
  quadtree.cpp)

enable_testing ()

if (BUILD_TESTS)
  add_subdirectory (tests)
endif (BUILD_TESTS)

install (TARGETS smartquadtree
  EXPORT quadtree
  DESTINATION lib
  INCLUDES DESTINATION include)

install (FILES neighbour.h quadtree.h quadtree.hpp
  DESTINATION include)

install (EXPORT quadtree
  DESTINATION share/cmake
  FILE quadtree-config.cmake)

