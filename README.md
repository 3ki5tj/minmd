# MinMD

A minimalist molecular dynamics (MD) simulation package.

## Overview

Our aim is to write a simple, didactic molecular dynamics program that helps people
understand the main logic behind many MD algorithms.
We wish to keep the code as simple as possible and the optimization is limited.

## Compiling 

Quick way of building the system:
```
make
```
The executables are under various binary directories,
we can run the MD simulation on the Lennard-Jones fluid by
```
./build/src/bin/md_lj/md_lj
```

For the debug version
```
make debug
valgrind ./debug/src/bin/md_lj/md_lj
```

See INSTALL.md for more details.


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


