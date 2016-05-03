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

/* ----------------------------------------------------------------------
   this code is used to find the neighbors of channel elements (square lattice)
------------------------------------------------------------------------- */

#include "string.h"
#include "math.h"
#include "deformation.h"
#include "group.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "stdlib.h" 
#include "math.h" 
#include "mpi.h"
#include "comm.h"
#include "update.h"
#include "neighbor.h"
#include "memory.h"
#include "domain.h"
#include "timer.h"
#define TINY  1.e-3 ;
#define FAKE_INT_VALUE  -991 ;

using namespace LAMMPS_NS;


double global_sxx;
double global_syy;
double global_szz;
double global_angle;
double global_young;
double global_poisson;

/* ---------------------------------------------------------------------- */

Deformation::Deformation(LAMMPS *lmp) : Pointers(lmp) {}

void Deformation::command(int narg, char **arg)
{

  if (narg < 4) error->all(FLERR,"Illegal deformation command"); 
  /*
  double  xfactor = atof(arg[0]);
  double  yfactor = atof(arg[1]);
  double  zfactor = atof(arg[2]);
  */
  
  sxx = atof(arg[0]);  
  syy = atof(arg[1]);
  szz = atof(arg[2]);
  angle = 3.14159*(atof(arg[3]))/180; // counterclockwise angle in radians. 
                                              
  global_sxx = sxx;
  global_syy = syy;
  global_szz = szz;
  global_angle = angle;
  
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int n ;
  
  int **bond_atom = atom -> bond_atom;
  int *tag = atom->tag;
  double xtemp, ytemp, ztemp;
  double xxx, yyy, zzz;
  
  double Y_mod = 50000000000;  
  double P_rat = 0.25;

  global_young = Y_mod;
  global_poisson = P_rat; 


  /*
  for (n = 0; n < nlocal; n++){
    if (tag[n] ==681){ fprintf(screen, "b==bond_rock1,rock2, rock3, rock4 %d %d %d %d\n",bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3]); }
  }
  */
  for (n = 0; n < nlocal; n++){
    /*
    x[n][0] = xfactor* x[n][0];
    x[n][1] = yfactor* x[n][1];
    x[n][2] = zfactor* x[n][2];
    */
        
    xtemp = x[n][0];
    ytemp = x[n][1];
    ztemp = x[n][2];
    //    fprintf(screen,"the %d-th atom. xfactor yfactor %f %f x y %f %f \n", n,xfactor,yfactor, x[n][0],x[n][1]);

    xxx = xtemp;     // For rotation calculation in x-y plane
    yyy = ytemp;
    //    zzz = ztemp;
    
    xtemp = cos(angle) * xxx + sin(angle) * yyy;
    ytemp = -1*sin(angle) * xxx + cos(angle) * yyy;

    xtemp = xtemp * (1 + (sxx - P_rat*syy - P_rat*szz)/Y_mod );
    ytemp = ytemp * (1 + (-P_rat*sxx + syy - P_rat*szz)/Y_mod );
    ztemp = ztemp * (1 + (-P_rat*sxx - P_rat*syy + szz)/Y_mod );

    xxx = xtemp;      // For rotation calculation in x-y plane
    yyy = ytemp;
    //    zzz = ztemp;    

    xtemp = cos(angle) * xxx - sin(angle) * yyy;
    ytemp = sin(angle) * xxx + cos(angle) * yyy;
    
    x[n][0] = xtemp;
    x[n][1] = ytemp;
    x[n][2] = ztemp;    

    

  }

  /*

  for (n = 0; n < nlocal; n++){
    if (tag[n] ==681){
      fprintf(screen, "d==bond_rock1,rock2, rock3, rock4 %d %d %d %d\n",bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3]);
    }
  }
  */
  timer->stamp();
  comm->forward_comm();
  timer->stamp(TIME_COMM);
}

