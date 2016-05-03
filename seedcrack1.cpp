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
#include "seedcrack1.h"
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

#define TINY  1.e-3 ;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Seedcrack1::Seedcrack1(LAMMPS *lmp) : Pointers(lmp) {}

void Seedcrack1::command(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Illegal injection_ini command"); // ID group-ID fixchannelini

  seedcrack_x = atof(arg[0]);
  seedcrack_y = atof(arg[1]);
  seedcrack_z = atof(arg[2]);
  //  ini_channel_length = atof(arg[3]);
  ini_channel_width = atof(arg[3]);


  int n;
  int nlocal = atom->nlocal;
  double **x1 = atom->x1;
  double **x = atom->x;


  for (n =0; n<nlocal; n++){
    x1[n][0] = x[n][0];
    x1[n][1] = x[n][1];
    x1[n][2] = x[n][2];

  }
  
  insert_seedcrack();
  comm->forward_comm();
}

/*--------------------------------------*/
void Seedcrack1::insert_seedcrack()
{
  int nlocal = atom->nlocal;
  double **x0 = atom->x0;
  int *tag = atom->tag;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int newton_bond = force->newton_bond;
  int **bond_type = atom->bond_type;
 
  int temp_bond_atom; 
  double atom1x, atom1y, atom1z, atom2x, atom2y,atom2z, channelx, channely,channelz;
  int atom1,atom2,i;
  int   nmax =0;
  int j;
  int nall = atom->nlocal + atom->nghost;

  int type;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int n;
  int *atype = atom->type;
  double *channel_width = atom->channel_width;
  int m;
  int rock_atom1, rock_atom2,channel_atom;
  

  for (n=0; n<nlocal;n++){
    if (atype[n] != CHANNEL_ATOM_TYPE) continue;
    channelx = x0[n][0];
    channely = x0[n][1];
    channelz = x0[n][2];
    if ((fabs(channelx - seedcrack_x) < 1.e-8 ) && ( fabs(channely-seedcrack_y)<1.e-8 ) && fabs(channelz-seedcrack_z)<1.e-8)
      {
	atype[n] = CONNECTED_CHANNEL_ATOM_TYPE;
	channel_width[n] = ini_channel_width;
      }
  }


  partner = NULL;

  if (atom->nmax > nmax) {
    memory->destroy(partner);
    nmax = atom->nmax;
    memory->create(partner,nmax,"fix_injection_update:partner");
  }

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
  }
  for (n = 0; n < nbondlist; n++) {
    atom1 = bondlist[n][0];
    atom2 = bondlist[n][1];
    type = bondlist[n][2];
    
    atom1x=x0[atom1][0];
    atom1y=x0[atom1][1];
    atom1z=x0[atom1][2];
    atom2x=x0[atom2][0];
    atom2y=x0[atom2][1];
    atom2z=x0[atom2][2];
    channelx=0.5*(atom1x+atom2x);
    channely=0.5*(atom1y+atom2y);
    channelz=0.5*(atom1z+atom2z);
    
    if ((fabs(channelx - seedcrack_x) < 1.e-8 ) && ( fabs(channely-seedcrack_y)<1.e-8 ) && ( fabs(channelz-seedcrack_z)<1.e-8 ))
      { 
	partner[atom1] = tag[atom2];
	partner[atom2] = tag[atom1];
      }
  }
    
  
  for (i = 0; i < nlocal; i++) {
    if (partner[i] == 0) continue;
    j = atom->map(partner[i]);
    //    if (partner[j] != tag[i])  {fprintf(screen,"!!! something wrong!!!\n"); continue;}
  
    for (m = 0; m < num_bond[i]; m++ ){

      if (bond_atom[i][m] == partner[i]){
	bond_type[i][m] = WET_CHANNEL_BOND_TYPE; //this is a wet channel
      }
    }
  }



}

/*------------------------------------------------------*/
int Seedcrack1::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
 int i,j,m;

  m = 0;
  
  
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = partner[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void Seedcrack1::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){
    partner[i] = static_cast<int> (buf[m++]);
  }
}
/* ---------------------------------------------------------------------- */

int Seedcrack1::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){
    buf[m++] = partner[i];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void Seedcrack1::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    partner[j] += static_cast<int> (buf[m++]);
  }
}


/*=============================

  =============================*/

double Seedcrack1::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*sizeof(int);
  return bytes;
}

/*=============================

  =============================*/
int Seedcrack1::find_channel_atom(int rock_atom1, int rock_atom2){
  
  int ii,jj;
  int return_channel_atom;
  double channelx,channely,channelz;
  double **x0 = atom->x0;
  int *num_bond_rock_channel_atom = atom->num_bond_rock_channel_atom;
  int **bond_rock_channel_atom = atom-> bond_rock_channel_atom;
  int *num_bond_channel_atom = atom->num_bond_channel_atom;
  int channel_atom;
  
  
  return_channel_atom = -991;
  for (ii = 0; ii<num_bond_rock_channel_atom[rock_atom1]; ii++){
    
    channelx = 0.5* ( x0[rock_atom1][0] + x0[rock_atom2][0]);
    channely = 0.5* ( x0[rock_atom1][1] + x0[rock_atom2][1]);
    channelz = 0.5* ( x0[rock_atom1][2] + x0[rock_atom2][2]);
    
    channel_atom = atom->map(bond_rock_channel_atom[rock_atom1][ii]);
    if ( (fabs(x0[channel_atom][0] - channelx) <1.e-6) &&  (fabs(x0[channel_atom][1] - channely) <1.e-6)  &&  (fabs(x0[channel_atom][2] - channelz) <1.e-6)  ){ 
      //	atype[channel_atom] = 3;
      return_channel_atom = channel_atom;
      //      fprintf(screen, "======num channel atom is %d channel_atom %d, channel_atomx, channel_atomy %f %f \n",num_bond_rock_channel_atom[channel_atom],channel_atom, x0[channel_atom][0],x0[channel_atom][1]);  
      //      fprintf(screen, "+++ num bond channel atom is %d \n",num_bond_channel_atom[channel_atom]);
    }
  }
  return atom->map(return_channel_atom);
}



