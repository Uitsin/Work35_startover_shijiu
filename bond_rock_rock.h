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

CommandStyle(bond_rock_rock,BondRockRock)

#else

#ifndef LMP_BOND_ROCK_ROCK_H
#define LMP_BOND_ROCK_ROCK_H

#include "pointers.h"

namespace LAMMPS_NS {

class BondRockRock : protected Pointers {
 public:

  BondRockRock(class LAMMPS *);
  void command(int, char **);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  int * sort( int *,int);


 private:

  int *partner0;
  int *partner1;
  int *partner2;
  int *partner3;
  int *partner4;
  int *partner5;


};

}

#endif
#endif
