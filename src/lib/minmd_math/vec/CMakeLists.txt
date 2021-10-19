add_library(vec INTERFACE)

add_subdirectory(vec2d)
add_subdirectory(vec3d)

get_property(DIM GLOBAL PROPERTY MINMD_DIM)
if (DIM EQUAL 2)
  set(real_vec vec2d)
else()
  set(real_vec vec3d)
endif()

target_link_libraries(vec INTERFACE ${real_vec})
target_include_directories(vec
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# test.c
add_executable(test_vec test.c)
target_link_libraries(test_vec PRIVATE m vec)

