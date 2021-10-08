add_library(vel INTERFACE)

target_include_directories(vel
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test program
add_executable(test_vel test.c)
target_link_libraries(test_vel m)
target_include_directories(test_vel
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
  )

