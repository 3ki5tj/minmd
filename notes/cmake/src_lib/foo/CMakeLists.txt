project(foo)

add_subdirectory(def) # a header-only child library
add_subdirectory(bar) # a child library that contains source files

# FOO_SOURCES contains .c or .cpp files
# this variable is set by the CMakefiles.txt in child directories
# child libraries are interface libraries
# so their source files are forwarded to the parent directory here
message("foo sources: ${FOO_SOURCES}") 

add_library(${PROJECT_NAME} ${FOO_SOURCES})

target_link_libraries(${PROJECT_NAME} bar)

target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)



