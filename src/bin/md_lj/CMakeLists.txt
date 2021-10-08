project(md_lj)
add_executable(${PROJECT_NAME} md_lj.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE libminmd)

