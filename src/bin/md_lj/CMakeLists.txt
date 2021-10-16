project(md_lj)
add_executable(${PROJECT_NAME} md_lj.c)
target_link_libraries(${PROJECT_NAME} PRIVATE minmd_lib)

