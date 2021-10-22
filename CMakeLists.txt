cmake_minimum_required(VERSION 3.6)

# set the project name and version
project(MinMD VERSION 0.1.0)

# activate this option by `cmake .. -DMINMD_DOUBLE=ON`
option(MINMD_DOUBLE "Use double precision floating-point math" OFF)
if (MINMD_DOUBLE)
  message("Using double-precision floating-point numbers")
else()
  message("Using single-precision floating-point numbers")
endif()

configure_file(MinMDConfig.h.in MinMDConfig.h)

add_subdirectory(src)

