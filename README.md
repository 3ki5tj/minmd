# MinMD

A minimalist molecular dynamics (MD) simulation package.

## Overview

The aim of this package is to provide a simple, didactic molecular dynamics program
that helps students understand the main logic behind various MD algorithms.

One of our design goals is to keep the code as simple as possible,
and keep technical complexities at bay.
To this end, we choose to
* Use C instead of C++ as the programming language.
* Limit the amount of optimizations

## Compiling

There are two ways of building the system.

1. CMake way: for building the all programs of the system
2. In-place way: for developers to build test programs for a single module

See INSTALL.md for more details.

### CMake build

Quick way of building the system:
```
make
```
This will build all executables, including test programs for all modules. 

The executables are under various binary directories.
For example, the MD simulation for the Lennard-Jones fluid is given by
```
./build/src/bin/md_lj/md_lj
```

#### Debug build

For the debug version
```
make debug
valgrind ./debug/src/bin/md_lj/md_lj
```

### In-place build

For test programs under the source code directory, there is
another quick-and-dirty way to build them.

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


## References

Various books and existing MD packages are helpful.

### MD Books
Understanding Molecular Simulation: From Algorithms to Applications (2nd Edition)
by Daan Frenkel and Berend Smit, Academic Press

### MD packages

GROMACS is probably the fastest CPU MD code on earth
https://github.com/gromacs/gromacs

But GROMACS has become very complex, the early version v1.6 has a much simpler code structure.
https://github.com/gromacs/gromacs/tree/release-1-6/src/gmxlib


