project(minmd_mdutils)

add_library(${PROJECT_NAME} INTERFACE)

# define the module's dependences
target_link_libraries(${PROJECT_NAME}
  INTERFACE
  def rng minmd_utils)

# make this module's headers available to other modules
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test program
add_executable(test_mdutils test.c)
target_link_libraries(test_mdutils ${PROJECT_NAME} -lm)

