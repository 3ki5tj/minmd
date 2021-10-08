project(libminmd)

add_subdirectory(def)
add_subdirectory(math)

# LIBMINMD_SOURCES would include all modules with .c or .cpp files
# header-only modules are excluded
set(LIBMINMD_SOURCES)
list(APPEND LIBMINMD_SOURCES
  ${MATH_SOURCES}
)
add_library(${PROJECT_NAME} ${LIBMINMD_SOURCES})


# Here, we include header-only modules and binary modules
# NOTE: the module "def" should be implicitly included in other modules, such as "vec" etc.
target_link_libraries(${PROJECT_NAME}
  math
  )



