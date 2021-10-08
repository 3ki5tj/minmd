add_library(def INTERFACE)

target_include_directories(def
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}> # root of build directory, where MinMDConfig.h resides
)

# test.c
add_executable(test_def test.c)
target_link_libraries(test_def def -lm)
