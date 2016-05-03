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
#include "delta_cutoff_randomness3.h"
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
#include "types.h"
#define TINY  1.e-3 ;
#define FAKE_INT_VALUE  -991 ;

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

DeltaCutoffRandomness3::DeltaCutoffRandomness3(LAMMPS *lmp) : Pointers(lmp) {}

void DeltaCutoffRandomness3::command(int narg, char **arg)
{

  if (narg < 6) error->all(FLERR,"Illegal delta_cutoff_randomness command"); 
  cutoff_mean_x = atof(arg[0]);
  cutoff_dev_x = atof(arg[1]);

  cutoff_mean_y = atof(arg[2]);
  cutoff_dev_y = atof(arg[3]);

  cutoff_mean_z = atof(arg[4]);
  cutoff_dev_z = atof(arg[5]);


  double random_value;
  /* initialize random seed: */
  srand(time(NULL));
  
  double *channel_delta_cutoff = atom -> channel_delta_cutoff;
  int nlocal = atom->nlocal;
  int n ,channel_atom;
  double **x0 = atom -> x0;
  double channel_x, channel_y, channel_z;
  int check_x, check_y, check_z;
  int *atype = atom->type;


  for (n = 0; n < nlocal; n++){
    if (atype[n] == ROCK_ATOM_TYPE) continue; // this is not a channel element
    channel_atom = n;
    channel_x = x0[channel_atom][0];
    channel_y = x0[channel_atom][1];
    channel_z = x0[channel_atom][2];
    
    check_x = (fmod(fabs(channel_x-0.5), 2.0) ==0); //and (flag_x == 1);
    check_y = (fmod(fabs(channel_y-0.5), 2.0) ==0); //and (flag_y == 1);
    check_z = (fmod(fabs(channel_z-0.5), 2.0) ==0); // and (flag_z == 1);
    
    channel_delta_cutoff[channel_atom] = 1000000; 	// make all bond unbreakable 

    /* generate random number between -0.5 and 0.5: */
    if (check_x){
      random_value =  (double)rand() / (RAND_MAX) -0.5 ;
      channel_delta_cutoff[n] = cutoff_mean_x + cutoff_dev_x * random_value;
    }
    if (check_y){
      random_value =  (double)rand() / (RAND_MAX) -0.5 ;
      channel_delta_cutoff[n] = cutoff_mean_y + cutoff_dev_y * random_value;
    }
    if (check_z){
      random_value =  (double)rand() / (RAND_MAX) -0.5 ;
      channel_delta_cutoff[n] = cutoff_mean_z + cutoff_dev_z * random_value;
    }
  }
  
  
 comm->forward_comm();


  //  timer->stamp();
  //  comm->forward_comm();
  //  timer->stamp(TIME_COMM);
}

