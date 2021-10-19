project(minmd_lib)

add_subdirectory(minmd_basic)
add_subdirectory(minmd_math)
add_subdirectory(minmd_md)

# LIBMINMD_SOURCES would include all modules with .c or .cpp files
# header-only modules are excluded
set(LIBMINMD_SOURCES)
list(APPEND LIBMINMD_SOURCES
  ${MATH_SOURCES}
)
add_library(${PROJECT_NAME} ${LIBMINMD_SOURCES})

target_link_libraries(${PROJECT_NAME}
  minmd_basic
  minmd_math
  minmd_md
  )

target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)



