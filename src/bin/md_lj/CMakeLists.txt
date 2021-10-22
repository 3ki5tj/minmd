project(md_lj)
add_executable(${PROJECT_NAME} md_lj.c)
target_link_libraries(${PROJECT_NAME} PRIVATE minmd_lib)

add_executable(${PROJECT_NAME}_etot md_lj_etot.c)
target_link_libraries(${PROJECT_NAME}_etot PRIVATE minmd_lib)

