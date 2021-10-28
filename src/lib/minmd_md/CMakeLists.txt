project(minmd_md)

add_subdirectory(mdutils)
add_subdirectory(thermostat)
add_subdirectory(ewald)

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME}
  INTERFACE 
  mdutils
  thermostat
  ewald)

target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

