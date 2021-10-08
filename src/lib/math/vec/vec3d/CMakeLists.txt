add_library(vec3d INTERFACE)
target_link_libraries(vec3d INTERFACE def)
target_include_directories(vec3d
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test.c
add_executable(test_vec3d test.c)
target_link_libraries(test_vec3d PRIVATE m def)
