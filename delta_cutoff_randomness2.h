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

#ifdef COMMAND_CLASS

CommandStyle(delta_cutoff_randomness2,DeltaCutoffRandomness2)

#else

#ifndef LMP_DELTA_CUTOFF_RANDOMNESS2_H
#define LMP_DELTA_CUTOFF_RANDOMNESS2_H

#include "pointers.h"

namespace LAMMPS_NS {

class DeltaCutoffRandomness2 : protected Pointers {
 public:

  DeltaCutoffRandomness2(class LAMMPS *);
  void command(int, char **);

  int find_channel_atom(int, int);


 private:
  double cutoff_mean, cutoff_dev;
  //  double injection_x, injection_y,injection_z;
  double factor_x, factor_y, factor_z;


};

}

#endif
#endif
