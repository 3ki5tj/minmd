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

# determine if we need a 2D or 3D system
# activate this option by `cmake .. -DMINMD_2D=ON`
option(MINMD_2D "Work with two-dimensional system" OFF)
if (MINMD_2D)
  set_property(GLOBAL PROPERTY MINMD_DIM 2)
else()
  set_property(GLOBAL PROPERTY MINMD_DIM 3)
endif()
get_property(DIM GLOBAL PROPERTY MINMD_DIM)
message("Working with ${DIM}-dimensional systems")

option(MINMD_MYMATH "Use homemade math functions" OFF)

configure_file(MinMDConfig.h.in MinMDConfig.h)

add_subdirectory(src)

