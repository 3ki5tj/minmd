#include "constraint.h"


constraint_t *init_settle(int water_model, int nmol,
    real (*x0)[DIM], real (*x1)[DIM], real (*v)[DIM])
{
  int *id_o, imol;

  XNEW(id_o, nmol);
  for (imol = 0; imol < nmol; imol++) {
    id_o[imol] = imol*3;
  }

  constraint_t *cs = constraint_settle_water_new_ez(water_model,
      nmol, id_o, x0, x1, v);

  free(id_o);

  return cs;
}


void print_vec(real *v)
{
  printf("    %+8.5f %+8.5f %8.5f\n", v[0], v[1], v[2]);
}


void test_stock(void)
{
  real x0[3][3] = {{0}},
       x1[3][3] = {{0}},
       v[3][3] = {{0}};

  constraint_t *cs = init_settle(CONSTRAINT_SETTLE_WATER_SPC,
      1, x0, x1, v);

  rng_t *rng = rng_new(RNG_TYPE_PCG, time(NULL));

  // initialize the water configuration
  constraint_settle_data_t *sed = (constraint_settle_data_t *) cs->data->idata;
  real dist_hn = sed->dist_hh*0.5;
  real dist_on = sed->dist_on;

  real x0_[3][DIM] = {
    {0, 0, 0},
    {-dist_on,  dist_hn, 0},
    {-dist_on, -dist_hn, 0},
  };

#define STOCK_NMOL 3

  real x1_stock[STOCK_NMOL][3][3] = {
    // symmetric model, the second rotation is not necessary 
    {{0.02, -0.0, 0.},
     {-dist_on+0.02,  dist_hn+0.05, 0.03},
     {-dist_on+0.02, -dist_hn-0.05, 0.03}},

    // a 90-degree rotation in the x-y plane
    {{0.0, 0.0, 0.},
     { dist_hn,  dist_on, 0.0},
     {-dist_hn,  dist_on, 0.0}},

    // generic example
    {{0.0, 0.0, 0.},
    {-dist_on+0.01,  dist_hn+0.02, 0.03},
    {-dist_on+0.07, -dist_hn+0.05, 0.06}},
  };
  
  int imol;
 
  for (imol = 0; imol < STOCK_NMOL; imol++) {
    vec_copy(x0[0], x0_[0]);
    vec_copy(x0[1], x0_[1]);
    vec_copy(x0[2], x0_[2]);

    vec_copy(x1[0], x1_stock[imol][0]);
    vec_copy(x1[1], x1_stock[imol][1]);
    vec_copy(x1[2], x1_stock[imol][2]);
    
    vec_rand_gauss(v[0], 0.2, rng);
    vec_rand_gauss(v[1], 1.0, rng);
    vec_rand_gauss(v[2], 1.0, rng);

    printf("BEFORE: %g %g %g\n", vec_dist(x1[0], x1[1]), vec_dist(x1[0], x1[2]), vec_dist(x1[1], x1[2]));
    print_vec(x1[0]);
    print_vec(x1[1]);
    print_vec(x1[2]);
    
    constraint_apply(cs, CONSTRAINT_APPLY_COORDINATES|CONSTRAINT_APPLY_VELOCITIES);

    printf("AFTER: %g %g %g\n", vec_dist(x1[0], x1[1]), vec_dist(x1[0], x1[2]), vec_dist(x1[1], x1[2]));
    print_vec(x1[0]);
    print_vec(x1[1]);
    print_vec(x1[2]);

    real dxpr[3][DIM], dvpr[3][DIM];
    printf("VELOCITY: %g %g %g\n",
        vec_dot(vec_diff(dxpr[2], x1[0], x1[1]), vec_diff(dvpr[2], v[0], v[1])),
        vec_dot(vec_diff(dxpr[0], x1[1], x1[2]), vec_diff(dvpr[0], v[1], v[2])),
        vec_dot(vec_diff(dxpr[1], x1[2], x1[0]), vec_diff(dvpr[1], v[2], v[0])));
    constraint_settle_print_vec_(v[0], "v0");
    constraint_settle_print_vec_(v[1], "v1");
    constraint_settle_print_vec_(v[2], "v2");

    printf("\n\n");
  }
  
  constraint_delete(cs);

  rng_delete(rng);
}


void init_water_config(constraint_t *cs,
    real (*x0)[DIM], real (*x1)[DIM], real (*v)[DIM])
{
  rng_t *rng = rng_new(RNG_TYPE_PCG, time(NULL));

  // initialize the water configuration
  constraint_settle_data_t *sed = (constraint_settle_data_t *) cs->data->idata;
  real dist_hn = sed->dist_hh*0.5;
  real dist_on = sed->dist_on;

  real x0_[3][DIM] = {
    {0, 0, 0},
    {-dist_on,  dist_hn, 0},
    {-dist_on, -dist_hn, 0},
  };
  real dxo[3], dxh1[3], dxh2[3];

  int imol;

  for (imol = 0; imol < sed->nmol; imol++) {
    vec_copy(x0[imol*3],   x0_[0]);
    vec_copy(x0[imol*3+1], x0_[1]);
    vec_copy(x0[imol*3+2], x0_[2]);

    vec_rand_gauss(dxo,  0.002, rng);
    vec_rand_gauss(dxh1, 0.010, rng);
    vec_rand_gauss(dxh2, 0.010, rng);
    vec_add(x1[imol*3],   x0[imol*3],   dxo);
    vec_add(x1[imol*3+1], x0[imol*3+1], dxh1);
    vec_add(x1[imol*3+2], x0[imol*3+2], dxh2);

    vec_rand_gauss(v[imol*3],   0.2, rng);
    vec_rand_gauss(v[imol*3+1], 1.0, rng);
    vec_rand_gauss(v[imol*3+2], 1.0, rng);
  }

  rng_delete(rng);
}


int test(void)
{
  int nmol = 1;
  real (*x0)[DIM];
  real (*x1)[DIM];
  real (*v)[DIM];

  XNEW(x0, 3*nmol);
  XNEW(x1, 3*nmol);
  XNEW(v,  3*nmol);

  constraint_t *cs = init_settle(CONSTRAINT_SETTLE_WATER_SPC,
      nmol, x0, x1, v);

  init_water_config(cs, x0, x1, v);

  printf("BEFORE: %g %g %g\n", vec_dist(x1[0], x1[1]), vec_dist(x1[0], x1[2]), vec_dist(x1[1], x1[2]));
  
  constraint_apply(cs, CONSTRAINT_APPLY_COORDINATES|CONSTRAINT_APPLY_VELOCITIES);

  printf("AFTER: %g %g %g\n", vec_dist(x1[0], x1[1]), vec_dist(x1[0], x1[2]), vec_dist(x1[1], x1[2]));

  constraint_delete(cs);

  free(x0);
  free(x1);
  free(v);

  return 0;
}

int main(void)
{
  test_stock();
  //test();
}

