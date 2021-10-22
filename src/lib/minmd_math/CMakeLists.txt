project(minmd_math)

add_subdirectory(rng)
add_subdirectory(vec)
add_subdirectory(stat_accum)

# INTERFACE means a library that has no real object
add_library(${PROJECT_NAME} INTERFACE)

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

