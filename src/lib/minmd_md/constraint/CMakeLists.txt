project(constraint)

add_library(${PROJECT_NAME} INTERFACE)

# define the module's dependences
target_link_libraries(${PROJECT_NAME}
  INTERFACE
  rng
  utils
  mdutils)

# make this module's headers available to other modules
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test programs
add_subdirectory(tests)

#add_executable(test_${PROJECT_NAME} test.c)
#target_link_libraries(test_${PROJECT_NAME} ${PROJECT_NAME} -lm)

