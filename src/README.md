# MinMD source code directory

## Overview

```
 +- lib                 libraries
 |  |
 |  +- minmd_basic      basic definitions
 |  |   |
 |  |   +- def          definitions
 |  |   +- utils        utility functions 
 |  |
 |  +- minmd_math       mathematical
 |  |   |
 |  |   +- rng          random number generator
 |  |   +- vec          3D/2D vector
 |  |   +- stat_accum   statistical accumulator
 |  |
 |  +- minmd_md         MD-related 
 |      |
 |      +- mdutils      MD-related utility 
 |      +- thermostat   thermostats
 |
 +- bin                 MD programs for different systems
    |
    +- md_lj            MD program for a Lennard-Jones fluid
```

libraries are all headers
