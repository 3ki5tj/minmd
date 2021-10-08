add_library(rng INTERFACE)
target_link_libraries(rng INTERFACE def)
target_include_directories(rng
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test.c
add_executable(test_rng test.c)
target_link_libraries(test_rng m rng)

