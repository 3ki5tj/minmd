project(def)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}> # root of build directory, where MinMDConfig.h resides
)

# test.c
add_executable(test_${PROJECT_NAME} test.c)
target_link_libraries(test_${PROJECT_NAME} ${PROJECT_NAME} -lm)
