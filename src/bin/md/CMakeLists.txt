project(md)
add_executable(${PROJECT_NAME} md.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE libminmd)

