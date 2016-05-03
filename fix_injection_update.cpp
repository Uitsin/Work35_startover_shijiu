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

#include "string.h"
#include "math.h"
#include "fix_injection_update.h"
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
#include "memory.h"
#include "neighbor.h"
#include "types.h"
#include "domain.h"
#include "lattice.h"
#define TINY  1.e-3 ;
using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixInjectionUpdate::FixInjectionUpdate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix channel/update command"); // ID group-ID fixchannelini
  
  injection_x = atof(arg[3]);
  injection_y = atof(arg[4]);
  injection_z = atof(arg[5]);
  injection_rate = atof(arg[6]);// bbl/min

/*
  int check_injection_point;
  int check_x, check_y, check_z;
  check_injection_point = 0;
  
  check_x = (fmod(fabs(injection_x-0.5), 2.0) ==0); //and (flag_x == 1);
  check_y = (fmod(fabs(injection_y-0.5), 2.0) ==0); //and (flag_y == 1);
  check_z = (fmod(fabs(injection_z-0.5), 2.0) ==0); // and (flag_z == 1);
  if (check_x == 1 && check_y == 0 && check_z == 0){check_injection_point = 1;} 
  else if (check_x == 0 && check_y == 1 && check_z == 0){check_injection_point = 1;} 
  else if (check_x == 0 && check_y == 0 && check_z == 1){check_injection_point = 1;} 
  else{
    fprintf(screen, "!! Wrong!!!The injection position is Wrong !!!!, make sure (x%2)=0.5 or (y%2)=0.5 or (z%2)=0.5 ");
    error->one(FLERR,"!! Wrong!!!The injection position is Wrong !!!!, make sure (x%2)=0.5 or (y%2)=0.5 or (z%2)=0.5 ");
  }
*/

  /*
  double H =30;//  meter
  injection_rate =  0.0088*(injection_rate/100.0)*(30.0/H); // m^3/sec
  */
  // one bbl/min = 1/377.39 m^3/s
  injection_rate =  1/377.39 *injection_rate; // m^3/sec
  

}

/*----------------------------------*/
FixInjectionUpdate::~FixInjectionUpdate()
{}


/*==========
  setmask
  ==========*/
int FixInjectionUpdate::setmask()
{
  int mask = 0;
   mask |= PRE_FORCE;
  return mask;
}


/*============================
  called before force routine
  ============================*/
void FixInjectionUpdate::pre_force(int vflag)
//  void FixInjectionUpdate::final_integrate()
{
  comm->forward_comm();
  injection_process();// injection water at injection position
  comm->forward_comm();

}



/*=================================================
  update new channel atom position after injection
  =================================================*/
void FixInjectionUpdate:: injection_process()
{
  int *atype = atom->type;
  int channel_atomi;
  double channelx,channely,channelz;
  int n;
  int nlocal = atom->nlocal;
  double **x0 = atom->x0;
  double *channel_width = atom->channel_width;
  int nstep = update->ntimestep;
  double dw;
  
  double lattice_mag = domain->lattice->xlattice;

  for (n=0; n<nlocal;n++){
    if (atype[n] != CONNECTED_CHANNEL_ATOM_TYPE) continue;
    channel_atomi = n;
    channelx = x0[channel_atomi][0];
    channely = x0[channel_atomi][1];
    channelz = x0[channel_atomi][2];
    if ((fabs(channelx - injection_x) > 1.e-8 ) || ( fabs(channely-injection_y)>1.e-8 ) ||  ( fabs(channelz-injection_z)>1.e-8 )) continue;// not injection position  

    dw = injection_rate * update->dt / (lattice_mag* lattice_mag);
    channel_width[channel_atomi] += dw;
    if ((update->ntimestep % 10000 )==0) fprintf(screen, "injection process: timestep %d, delta_w = %f (microns),channel_width %f (mm),  at x y z %f %f %f\n ",nstep,dw*1000000,channel_width[channel_atomi]*1000,x0[channel_atomi][0],x0[channel_atomi][1],x0[channel_atomi][2] );
  }
}

