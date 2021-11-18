# MinMD source code directory

## Overview

```
 +- lib                 libraries
 |  |
 |  +- minmd_basic      basic definitions and functions
 |  |   |
 |  |   +- def          definitions
 |  |   +- utils        utility functions 
 |  |
 |  +- minmd_math       mathematical functions
 |  |   |
 |  |   +- rng          random number generator
 |  |   +- vec          3D/2D vector
 |  |   +- stat_accum   statistical accumulator
 |  |
 |  +- minmd_md         MD-related 
 |      |
 |      +- mdutils      utilities 
 |      +- thermostat   thermostats
 |      +- ewald        Ewald sum
 |      +- constraint   constraint solvers
 |
 +- bin                 MD programs for different systems
    |
    +- md_lj            Lennard-Jones fluid
    +- md_ocp           One-component plasma
```

libraries are all headers
