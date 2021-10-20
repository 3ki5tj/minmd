project(bar)

# INTERFACE means a library that has no real object
add_library(${PROJECT_NAME} INTERFACE)

# collect source code files and define a variable BAR_SOURCES
# that is available to the parent directory
file(GLOB BAR_SOURCES *.c *.cpp)
list(FILTER BAR_SOURCES EXCLUDE REGEX test.*) # requires CMake 3.6
set(FOO_SOURCES ${BAR_SOURCES} PARENT_SCOPE)
message("\n\nFOO_SOURCES ${MATH_SOURCES}\n\n")

# files that this library depends
target_link_libraries(${PROJECT_NAME}
  INTERFACE 
  def
  -lm)

# make the headers available to users of this library
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

# test.cpp
add_executable(test test.cpp ${BAR_SOURCES})
target_link_libraries(test PRIVATE ${PROJECT_NAME})
