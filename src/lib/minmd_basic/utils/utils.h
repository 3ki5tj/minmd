#ifndef UTIL_H__
#define UTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "def.h"

#ifndef XNEW
#define XNEW(x, n) x = xnew_(n, sizeof(*(x)), #x)
INLINE void *xnew_(size_t n, size_t size, const char *name)
{
  void *p;

  if ((p = calloc(n, size)) == NULL) {
    fprintf(stderr, "no memory for variable '%s' %lu x %lu\n",
        name, (unsigned long) n, (unsigned long) size);
    exit(1);
  }
  return p;
}
#endif /* defined(XNEW) */



#ifndef XCLONE
#define XCLONE(x, y, size) x = xclone_(y, size, #x)
INLINE void *xclone_(void *y, size_t size, const char *name_x)
{
  void *x;

  x = xnew_(1, size, name_x);
  memcpy(x, y, size);

  return x;
}
#endif /* defined(XCLONE) */



#endif /* UTIL_H__ */

