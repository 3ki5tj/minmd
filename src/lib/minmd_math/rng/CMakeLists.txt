project(rng)

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE def)
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test.c
add_executable(test_${PROJECT_NAME} test.c)
target_link_libraries(test_${PROJECT_NAME} ${PROJECT_NAME} -lm)

