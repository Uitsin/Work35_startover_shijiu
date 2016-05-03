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
#include "delta_cutoff_randomness2.h"
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

DeltaCutoffRandomness2::DeltaCutoffRandomness2(LAMMPS *lmp) : Pointers(lmp) {}

void DeltaCutoffRandomness2::command(int narg, char **arg)
{

  if (narg < 5) error->all(FLERR,"Illegal delta_cutoff_randomness command"); 
  cutoff_mean = atof(arg[0]);
  cutoff_dev = atof(arg[1]);
  factor_x = atof(arg[2]);
  factor_y = atof(arg[3]);
  factor_z = atof(arg[4]);

  //  injection_x = atof(arg[2]);
  //  injection_y = atof(arg[3]);
  //  injection_z = atof(arg[4]);



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
      channel_delta_cutoff[n] = cutoff_mean + factor_x * cutoff_dev * random_value;
    }
    if (check_y){
      random_value =  (double)rand() / (RAND_MAX) -0.5 ;
      channel_delta_cutoff[n] = cutoff_mean + factor_y * cutoff_dev * random_value;
    }
    if (check_z){
      random_value =  (double)rand() / (RAND_MAX) -0.5 ;
      channel_delta_cutoff[n] = cutoff_mean + factor_z * cutoff_dev * random_value;
    }

    /*
    if (check_x || check_y || check_z){
      channel_delta_cutoff[n] = cutoff_mean;
    }

    if ((check_x && flag_x) || (check_y && flag_y) || (check_z && flag_z) ){
      //    if (check_x || check_y || check_z){
      random_value =  (double)rand() / (RAND_MAX) -0.5 ;
      //    fprintf(screen, "random value is %f \n ",random_value);    
      channel_delta_cutoff[n] = cutoff_mean + cutoff_dev * random_value;
    }
    else {
      fprintf(screen, "channel_delta_cutoff is  %f at x y z  %f %f %f and check_x, check_y, check_z %d %d %d \n  ", channel_delta_cutoff[channel_atom],channel_x, channel_y, channel_z,check_x, check_y, check_z);
    }
    */
    
    
    
  }
  
  

  comm->forward_comm();

  /*=============
    check if randomness routine works
    =============*/

  /*

  int m;
  int *atype = atom->type;
  int rock_atom1, rock_atom2, channel_atom;
  int *num_bond = atom-> num_bond;
  int **bond_atom = atom-> bond_atom;

  for (n =0; n <nlocal; n++) {
    if (atype[n] != 1) continue; // this is not a rock atom
    rock_atom1 = n;
    for (m = 0; m < num_bond[rock_atom1]; m++) {
      rock_atom2 = atom->map(bond_atom[rock_atom1][m]);
      if (atype[rock_atom2] !=1) {fprintf(screen, "warning!! this should be a rock atom!!\n"); continue;}
      channel_atom = find_channel_atom(rock_atom1, rock_atom2);
      fprintf(screen,"channel delta cutoff is %f\n", channel_delta_cutoff[channel_atom]);
    }
  }
  */
  //  timer->stamp();
  //  comm->forward_comm();
  //  timer->stamp(TIME_COMM);
}

/*====================================================
  find channel_atom between rock_atom1, and rock_atom2
  ====================================================*/

int DeltaCutoffRandomness2::find_channel_atom(int rock_atom1, int rock_atom2){
  
  int ii,jj;
  int return_channel_atom;
  double channelx,channely,channelz;
  double **x0 = atom->x0;
  int *num_bond_rock_channel_atom = atom->num_bond_rock_channel_atom;
  int **bond_rock_channel_atom = atom-> bond_rock_channel_atom;
  int *num_bond_channel_atom = atom->num_bond_channel_atom;
  int channel_atom;
  int *tag = atom->tag;
  
  
  return_channel_atom = -991;
  for (ii = 0; ii<num_bond_rock_channel_atom[rock_atom1]; ii++){
    
    channelx = 0.5* ( x0[rock_atom1][0] + x0[rock_atom2][0]);
    channely = 0.5* ( x0[rock_atom1][1] + x0[rock_atom2][1]);
    channelz = 0.5* ( x0[rock_atom1][2] + x0[rock_atom2][2]);
    
    channel_atom = atom->map(bond_rock_channel_atom[rock_atom1][ii]);
    if ( (fabs(x0[channel_atom][0] - channelx) <1.e-6) &&  (fabs(x0[channel_atom][1] - channely) <1.e-6) &&  (fabs(x0[channel_atom][2] - channelz) <1.e-6) ){ 
      return_channel_atom = channel_atom;
    }
  }
  return return_channel_atom;
}



