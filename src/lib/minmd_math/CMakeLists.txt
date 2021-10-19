project(minmd_math)

add_subdirectory(rng)
add_subdirectory(vec)
add_subdirectory(stat_accum)

# INTERFACE means a library that has no real object
add_library(${PROJECT_NAME} INTERFACE)

# collect source code files and define a variable MATH_SOURCES
# that is available to the parent directory
file(GLOB MATH_SOURCES *.c *.cpp)
#list(REMOVE_ITEM MATH_SOURCES
#    ${CMAKE_CURRENT_SOURCE_DIR}/test.cpp)
list(FILTER MATH_SOURCES EXCLUDE REGEX test.*) # requires CMake 3.6
set(MATH_SOURCES ${MATH_SOURCES} PARENT_SCOPE)
message("\n\nMATH_SOURCES ${MATH_SOURCES}\n\n")


target_link_libraries(${PROJECT_NAME}
  INTERFACE 
  rng
  vec
  stat_accum
  -lm)

# make the headers available to users of this library
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

# for building test_math
add_executable(test_${PROJECT_NAME} test.cpp ${MATH_SOURCES})
target_link_libraries(test_${PROJECT_NAME} PRIVATE ${PROJECT_NAME})
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
