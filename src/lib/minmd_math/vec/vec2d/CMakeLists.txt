add_library(vec2d INTERFACE)
target_link_libraries(vec2d INTERFACE def)
target_include_directories(vec2d
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test.c
add_executable(test_vec2d test.c)
target_link_libraries(test_vec2d PRIVATE m def)
