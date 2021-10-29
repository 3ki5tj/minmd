project(md_ocp)
add_executable(${PROJECT_NAME} md_ocp.c)
target_link_libraries(${PROJECT_NAME} PRIVATE minmd_lib)

