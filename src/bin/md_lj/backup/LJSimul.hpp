#include <iostream>
#include <vector>
#include <ctime>
#include "math.h"


class LJPotentialCutoff
{
  private:
    real rc_def_; // preferred cutoff distance
    real rc_; // actual cutoff distance

  public:
    LJPontentialCutoff(real rc_def, real l) :
      rc_def_(rc_def),
      rc_(rc_def < l*0.5 ? rc_def : l*0.5) {}

    real get(void) { return rc_; }
};


class LJPotentialEnergyData
{
  private:
    int dim_;
    int n_;
    real rc_;
    real rho_;
    real epot_ = 0,
         ep0_ = 0,
         ep6_ = 0,
         ep12_ = 0,
         eps_ = 0;
    real epot_shift_;
    real epot_tail_;
    real vir_ = 0;
    real p_tail_;

  public:
    LJPotentialEnergyData(dim, n, rc, rho) :
      dim_(dim),
      n_(n),
      rc_(rc),
      rho_(rho),
    {
      real irc = 1 / rc_,
           irc2 = irc * irc,
           irc4 = irc2 * irc2,
           irc6 = irc4 * irc2;

      epot_shift_ = 4 * irc6 * (irc6 - 1);
      if (dim_ == 2) {
        epot_tail_ = PI * rho_ * n_ * (0.4*irc6 - 1) * irc4;
        p_tail_ = PI * rho_ * rho_ * (2.4*irc6 - 3) * irc4;
      } else {
        epot_tail = 8 * PI * rho_ * n_ / 9 * (irc6 - 3) * irc3;
        p_tail = 32 * PI * rho_ * rho_ / 9 * (irc6 - 1.5) * irc3;
      }
    }
};


class LJSimul
{
  private:
    int dim_ = DIM; // dimension, 3 or 2
    int n_; // number of particles
    int n_dof_; // number of degrees of freedom
    real rho_;
    real vol_;
    real l_;
    LJPotentialCutOff cutoff_;

    vector<real [DIM]> x_; // position
    vector<real [DIM]> v_; // velocity
    vector<real [DIM]> f_; // force

    vector<vector<real>> r2ij_; // cached pair distances 

    LJPontetialEnergyData epot_data_;
    real ekin_; // kinetic energy

  public:
    LJSimul(int n, real rho, real rc_def);
    void do_md(long long int nsteps);
    void save_checkpoint(const char *fn);
    void print_summary(void);
};



LJSimul::LJSimul(int n, real rho, real rc_def)
  : dim_(DIM),
    n_(n),
    n_dof_(n*DIM - DIM),
    rho_(rho),
    vol_(n/rho),
    l_(pow(n/rho, 1.0/DIM)),
    cutoff_(rc_def, l_),
    x_(n),
    v_(n),
    f_(n),
    r2ij_(n, vector<real>(n)),
    epot_data_(DIM, rc_),
    ekin_(0)
{
  // initialize positions on the lattice
  // lj_init_fcc(lj);

  // initialize random velocities
  // v[i][d] = random_gauss()

  // remove center of mass motion
}



void LJSimul::do_md(long long int nsteps)
{
}

void LJSimul::save_checkpoint(const char *fn)
{
}

void LJSimul::print_summary(void)
{
}

