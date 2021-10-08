project(math)

add_subdirectory(rng)
add_subdirectory(vec)

# INTERFACE means a library that has no real object
add_library(${PROJECT_NAME} INTERFACE)

# collect source code files and define a variable MATH_SOURCES
# that is available to the parent directory
file(GLOB MATH_SOURCES *.cpp)
list(REMOVE_ITEM MATH_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/test.cpp)
set(MATH_SOURCES ${MATH_SOURCES} PARENT_SCOPE)

# this project depends on "def"
target_link_libraries(${PROJECT_NAME}
  INTERFACE 
  def rng vec)

# make the headers available to users of this library
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

# for building test_math
add_executable(test_math test.cpp ${MATH_SOURCES})
target_link_libraries(test_math PRIVATE def)
#
# Since `math` is an interface library,
# the library contains no binary objects
# and the .cpp files are forwarded to the parent level via ${MATH_SOURCES}
# so we cannot do something like this when building the executable
#
#   add_executable(test_math test.cpp)
#   target_link_libraries(test_math PRIVATE ${PROJECT_NAME})
#
# functions in the .cpp will be missing in the linking step
#
