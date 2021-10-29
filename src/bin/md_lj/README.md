## Lennard-Jones fluid

A box of particles interacting with one another via the 12-6 Lennard-Jones potential
$$
  u(r_{ij}) = 4 \left[ \left( \frac{\sigma}{r_{ij}) \right)^12 - \left( \frac{\sigma}{r_{ij}) \right)^6 \right].
$$

## Programs

* md_lj       main MD program
* md_lj_etot  test the conservation of total energy
* force_test  test the consistency of forces

## Building

Under the root directory of MinMD
```
make
```

Run the program as
```
./build/src/bin/md_lj/md_lj
```

## Testing

```
make md_lj && valgrind ./md_lj
```

NOTE: the Makfile under this directory is independent of the CMakeLists.txt (and hence the CMake system).
It is used to do compile the test programs.

