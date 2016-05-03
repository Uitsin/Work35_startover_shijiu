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

#ifdef FIX_CLASS

FixStyle(lubrication/update_BC,FixLubricationUpdateBC)

#else

#ifndef LMP_FIX_LUBRICATION_UPDATE_BC_H
#define LMP_FIX_LUBRICATION_UPDATE_BC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLubricationUpdateBC : public Fix {
 public:
  FixLubricationUpdateBC(class LAMMPS *, int, char **);
  ~FixLubricationUpdateBC();

  int setmask();

  //void pre_force(int);
  void final_integrate();

  int find_channel_atom(int, int);
  void cal_channel_pressure();
  void lubrication();
  void copy_channel_width();
  void check_channel_pressure();
  void channel_update();
  void new_channel();  
  void bond_break();

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  double BC_xlo, BC_xhi;
  double BC_ylo, BC_yhi;
  double BC_zlo, BC_zhi;
  double vij_max;
  //  double injection_x, injection_y, injection_z;
  int nmax;
  int *partner;

};

}

#endif
#endif
