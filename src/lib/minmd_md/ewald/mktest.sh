#!/usr/bin/bash
# for manually compiling test.c
gcc -g test.c -o test \
  -I../../../../debug \
  -I../../minmd_basic/def \
  -I../../minmd_basic/utils \
  -I../../minmd_math/rng \
  -I../../minmd_math/vec \
  -lm

