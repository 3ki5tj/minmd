project(minmd_lib)

add_subdirectory(def)
add_subdirectory(utils) # minmd_utils
add_subdirectory(math) # minmd_math
add_subdirectory(mdutils) # minmd_mdutils

# LIBMINMD_SOURCES would include all modules with .c or .cpp files
# header-only modules are excluded
set(LIBMINMD_SOURCES)
list(APPEND LIBMINMD_SOURCES
  ${MATH_SOURCES}
)
add_library(${PROJECT_NAME} ${LIBMINMD_SOURCES})

target_link_libraries(${PROJECT_NAME}
  minmd_utils
  minmd_math
  minmd_mdutils
  -lm
  )

target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)



