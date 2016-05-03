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

CommandStyle(seedcrack1,Seedcrack1)

#else

#ifndef LMP_SEEDCRACK1_H
#define LMP_SEEDCRACK1_H

#include "pointers.h"

namespace LAMMPS_NS {

class Seedcrack1 : protected Pointers {
 public:
  Seedcrack1(class LAMMPS *);

  void command(int, char **);
  void insert_seedcrack();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  int find_channel_atom(int,int);


 private:
  int **channel_neighbor;
  int *partner;
  double seedcrack_x, seedcrack_y, seedcrack_z;
  double ini_channel_length;
  double ini_channel_width;

};

}

#endif
#endif
