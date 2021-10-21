## Source code library example

In the `minmd_lib` library, we follow the policy of using header only files.
For reference, we give in `src_lib` an example of usual library using .c and .cpp files
See the README.md there for details.

## Common ommands

### Print internal variables

```
message("extra-libraries: ${EXTRA_LIBS}")
```

### set

https://cmake.org/cmake/help/latest/command/set.html

```
set(<variable> <value>... [PARENT_SCOPE])
```
Each new directory or function creates a new scope.

```
set(ENV{<variable>} [<value>])
```

### set_property

```
set_property(GLOBAL PROPERTY global_var_name 3)
```

retrieve to a local variable
```
get_property(local_var GLOBAL PROPERTY global_var_name)
```



### Difference between public, private and interface links

https://stackoverflow.com/questions/26037954/cmake-target-link-libraries-interface-dependencies

https://stackoverflow.com/a/67385469/3326606

When A links B as PRIVATE,
it is saying that A uses B in its implementation,
but B is not used in any part of A's public API.
Any code that makes calls into A would not need to refer directly to anything from B.
An example of this could be a networking library A which can be built to use
one of a number of different SSL libraries internally (which B represents).
A presents a unified interface for client code
which does not reference any of the internal SSL data structures or functions.
Client code would have no idea what SSL implementation (B) is being used by A,
nor does that client code need to care.

When A links B as INTERFACE, it is saying that
A does not use B in its implementation,
but B is used in A's public API.
Code that calls into A may need to refer to things from B in order to make such calls.
One example of this is an interface library which simply forwards calls along
to another library but doesn't actually reference the objects on the way
through other than by a pointer or reference.
Another example is where A is defined in CMake as an interface library,
meaning it has no actual implementation itself,
it is effectively just a collection of other libraries
(I'm probably over-simplifying here, but you get the picture).


INTERFACE means things that consumers require but the producer doesn't.
https://cmake.org/cmake/help/latest/guide/tutorial/Adding%20Usage%20Requirements%20for%20a%20Library.html


### Header-only library

```
+- header1
|   |
|   +- header1.h
|
+- header2
|   |
|   +- header2.h (uses header1.h)
|
+- app
    |
    +- app.c (uses header2.h) 
```

header1/CMakeLists.txt:
```
add_library(header1 INTERFACE)
target_include_directories(header1
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )
```

header2/CMakeLists.txt:
```
add_library(header2 INTERFACE)
target_link_libraries(header2 INTERFACE header1)
target_include_directories(header2
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )
```

app/CMakeLists.txt:
```
add_executable(app app.c)
target_link_libraries(app header2)
```


