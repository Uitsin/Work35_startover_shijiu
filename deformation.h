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

CommandStyle(deformation,Deformation)

#else

#ifndef LMP_DEFORMATON_H
#define LMP_DEFORMATON_H

#include "pointers.h"

extern double global_sxx;
extern double global_syy;
extern double global_szz;
extern double global_angle;
extern double global_young;
extern double global_poisson;

namespace LAMMPS_NS {

class Deformation : protected Pointers {
 public:

  Deformation(class LAMMPS *);
  void command(int, char **);

  double  sxx;  
  double  syy;
  double  szz;
  double  syz;
  double  sxz;
  double  sxy;
  double  angle; // counterclockwise angle in radians. 
  double  mix_start;


 private:



};

}

#endif
#endif
