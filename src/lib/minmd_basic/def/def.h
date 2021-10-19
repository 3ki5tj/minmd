#ifndef DEF_H__
#define DEF_H__

#include "MinMDConfig.h"

#ifndef INLINE
#define INLINE static inline
#endif

#ifdef MINMD_DOUBLE
typedef double real;
#else
typedef float real;
#endif

#ifdef MINMD_2D
#define DIM 2
#else
#define DIM 3
#endif

const real PI = 3.141592653589793;

#endif /* DEF_H__ */

