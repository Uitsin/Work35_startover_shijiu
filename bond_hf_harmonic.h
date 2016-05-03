/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(hf_harmonic,BondHfHarmonic)

#else

#ifndef LMP_HF_BOND_HARMONIC_H
#define LMP_HF_BOND_HARMONIC_H

#include "stdio.h"
#include "bond.h"

namespace LAMMPS_NS {

class BondHfHarmonic : public Bond {
 public:
  BondHfHarmonic(class LAMMPS *);
  virtual ~BondHfHarmonic();
  virtual void compute(int, int);
  void coeff(int,  char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int, double &);

  void rock_force_update();
  //  void channel_force_update();
  

 protected:
  double *r0, *rg, *k1,*k2, *beta;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
