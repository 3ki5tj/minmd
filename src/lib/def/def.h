#ifndef DEF_H__
#define DEF_H__

#include "MinMDConfig.h"

#ifdef MINMD_DOUBLE
typedef double real;
#else
typedef float real;
#endif

#ifdef MINMD_2D
#define MINMD_DIM 2
#else
#define MINMD_DIM 3
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef INLINE
#define INLINE static inline
#endif

#endif /* DEF_H__ */

