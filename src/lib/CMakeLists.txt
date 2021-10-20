project(minmd_lib)

add_subdirectory(minmd_basic)
add_subdirectory(minmd_math)
add_subdirectory(minmd_md)

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(${PROJECT_NAME}
  INTERFACE
  minmd_basic
  minmd_math
  minmd_md
  )

target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)



