## One-component plasma

This system is a box of charged Lennard-Jones fluids.
All particles carry the same charge.
And the system is immersed in a uniformly compenstating charge background.

Computationally, we only need to add Ewald sum to the system,
ignoring the k = 0 component.


## Building

Under the root directory of MinMD
```
make
```

Run the program as
```
./build/src/bin/md_ocp/md_ocp
```

## Testing

```
make md_ocp && valgrind ./md_ocp
```

NOTE: the Makfile under this directory is independent of the CMakeLists.txt (and hence the CMake system).
It is used to do compile the test programs.

