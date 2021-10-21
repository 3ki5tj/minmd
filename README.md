# MinMD

A minimalist molecular dynamics (MD) simulation package.

## Overview

Our aim is to write a simple, didactic molecular dynamics program that helps people
understand the main logic behind various MD algorithms.
We wish to keep the code as simple as possible and the optimization is limited.

## Compiling 

Quick way of building the system:
```
make
```
The executables are under various binary directories, e.g. `build/src/bin/md_lj/md_lj`

```
./build/src/bin/md_lj/md_lj
```

For the debug version
```
make debug
valgrind ./debug/src/bin/md_lj/md_lj
```

See INSTALL for more details.


## References

Various books and existing MD packages are helpful.

Understanding Molecular Simulation: From Algorithms to Applications (2nd Edition)
by Daan Frenkel and Berend Smit, Academic Press

GROMACS v1.6
https://github.com/gromacs/gromacs/tree/release-1-6/src/gmxlib

GROMACS v1.6 has a simple code base.
