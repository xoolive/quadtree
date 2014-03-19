cmake_minimum_required(VERSION 2.6)
project (detection CXX)

find_program (INTEL_COMPILER icc)
if (INTEL_COMPILER)
  include (CMakeForceCompiler)
  cmake_force_c_compiler (icc "Intel C Compiler")
  cmake_force_cxx_compiler (icpc "Intel C++ Compiler")
endif (INTEL_COMPILER)


enable_testing ()

add_executable (test_neighbour
  test_neighbour.cpp)

add_test (neighbour
  test_neighbour)

add_executable (test_quadtree
  test_quadtree.cpp
  quadtree.cpp)

add_test (quadtree
  test_quadtree)

find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)

option (RUN_IN_VM "Run inside a virtual machine")

# no-deprecated-declarations for OpenGL/GLU and GLUT
set (CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations")

include_directories (
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  )

add_executable (test_simu
  test_simu.cpp
  quadtree.cpp)

target_link_libraries(test_simu
  ${OPENGL_gl_LIBRARY}
  ${OPENGL_glu_LIBRARY}
  ${GLUT_glut_LIBRARY}
  )
 
file (COPY
  "run.cmake"
  DESTINATION ${CMAKE_BINARY_DIR}
  )

if (RUN_IN_VM)
  add_custom_target (run
    COMMAND ${CMAKE_COMMAND} -P run.cmake
    DEPENDS test_simu
    )
else(RUN_IN_VM)
  add_custom_target (run
    COMMAND test_simu
    DEPENDS test_simu
    )
endif(RUN_IN_VM)
