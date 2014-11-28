find_package (OpenGL REQUIRED)
find_package (GLUT REQUIRED)

option (RUN_IN_VM "Run inside a virtual machine")

if (APPLE)
  # no-deprecated-declarations for OpenGL/GLU and GLUT
  set (CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations")
endif (APPLE)

macro (prepare_test target)
  add_executable (test_${target} test_${target}.cpp)
  target_link_libraries (test_${target} smartquadtree)
  add_test (${target} test_${target})
endmacro (prepare_test)

prepare_test (neighbour)
prepare_test (quadtree)
prepare_test (clipping)

include_directories (
  ".."
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  )

add_executable (test_simu
  test_simu.cpp)

target_link_libraries(test_simu
  smartquadtree
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
