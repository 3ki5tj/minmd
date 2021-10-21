## Overview

We use cmake as our building system

## CMake

### Installing cmake

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
