## Overview

There are two ways of compiling the programs in this package

The first way, which is the official way, is to use CMake.
This will build all executables in the packages in one settings,
and the executables are saved under a separate directory from the source code tree.

The second way is the in-place build.
This is for developers to do quick testings for some individual module.
It relies on the Makefile under the each source code directory
to produce only test programs under that directory.

For example, if we type `make` under `src/lib/minmd_md/thermostat`,
it will compile `test.c` to a binary program `test` under that directory.

Note that under each directory, the CMakeLists.txt and Makefile are independent.
The CMakeLists.txt is used to used for the CMake build,
while the Makefile is used for making the in-place test program(s). 

## CMake build

### Installing CMake

On Ubuntu, we can install cmake by typing
```
sudo apt install cmake
```


## Building the project

### Automatical way

Release version:
```
make
```

Debug version:
```
make debug
```

### Manually building the project

Under the project root directory, where this file is located

```
mkdir build
cd build
cmake ..
cmake --build .
```

After the binaries are successfully build, the `md` program is located under `build/src/programs/md`

### CMake tutorials

https://cmake.org/cmake/help/latest/guide/tutorial/index.html

https://riptutorial.com/cmake

See also notes/cmake.md for some notes.


## In-place build

For test programs under the source code directory, there is
a quick way to compile them.

For example, to build the test program under `src/lib/minmd_md/thermostat`,
simply type
```
make
```

This will produce the test program `test` under the `thermostat` directory.
To run it directly
```
./test
```

Or check for memory leaks
```
valgrind ./test
```



