project(minmd_md)

add_subdirectory(mdutils)
add_subdirectory(thermostats)

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME}
  INTERFACE 
  mdutils
  thermostats)
target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

