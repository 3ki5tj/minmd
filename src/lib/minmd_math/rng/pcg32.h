#ifndef PCG_H__
#define PCG_H__

/* Permuted Congruential Generator
 * https://www.pcg-random.org/download.html
 * https://en.wikipedia.org/wiki/Permuted_congruential_generator */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

typedef struct {
  uint64_t state;
  uint64_t inc;
} pcg32_random_t;

pcg32_random_t *pcg32_random_new(uint64_t seed, uint64_t inc)
{
  pcg32_random_t *r = (pcg32_random_t *) malloc(sizeof(pcg32_random_t));
  if (r == NULL) {
    exit(1);
  }

  if (inc == 0) {
    inc = 1442695040888963407ULL;
  }
  inc |= 1; /* ensure the increment is odd */

  if (seed == 0) {
    seed = 0x4d595df4d0f33173ULL;
  } else {
    seed += inc;
  }

  r->state = seed;
  r->inc = inc;
  return r;
}


void pcg32_random_delete(pcg32_random_t *r)
{
  free(r);
}


uint32_t pcg32_random_uint32(pcg32_random_t* pcg)
{
  uint64_t oldstate = pcg->state;
  // advance internal state;
  pcg->state = oldstate * 6364136223846793005ULL + pcg->inc;
  // calculate output function (XSH RR), uses old state for max ILP
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}


int pcg32_random_save(pcg32_random_t *pcg, const char *fn)
{
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "PCG32 %" PRIu64 " %" PRIu64 "\n",
      pcg->state, pcg->inc); 
  fclose(fp);
  return 0;
}

int pcg32_random_load(pcg32_random_t *pcg, const char *fn)
{
  FILE *fp;
  char s[1024];

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  if (fgets(s, sizeof s, fp) != NULL) {
    uint64_t state = 0, inc = 1;
    if (sscanf(s, "PCG32 %" SCNu64 " %" SCNu64, &state, &inc) == 2) {
      pcg->state = state;
      pcg->inc = inc;
    }
  }
  fclose(fp);
  return 0;
}

#endif /* PCG_H__ */

