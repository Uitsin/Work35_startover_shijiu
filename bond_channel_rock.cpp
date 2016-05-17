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
#include "bond_channel_rock.h"
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
#include "lattice.h"
#define TINY  1.e-3 ;
#define FAKE_INT_VALUE  -991 ;

#include "types.h"
//const int ROCK_ATOM_TYPE = 1;
//const int CHANNEL_ATOM_TYPE = 2;


using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

BondChannelRock::BondChannelRock(LAMMPS *lmp) : Pointers(lmp) {}

void BondChannelRock::command(int narg, char **arg)
{

  //  if (narg < 3) error->all(FLERR,"Illegal injection_ini command"); // ID group-ID fixchannelini
 

  int nmax =0;
  int nall = atom->nlocal + atom->nghost;

  int n, m,i;
  int *atype = atom->type;
  int *tag = atom->tag;
  int atomi, atomj;
  double xi,yi,zi;
  double xj,yj,zj;
  double **x0 = atom->x0;
  int nlocal = atom->nlocal;
  int cond0,cond1,cond2,cond3,cond4,cond5;
  int **bond_rock_channel_atom = atom->bond_rock_channel_atom;
  int *num_bond_rock_channel_atom = atom->num_bond_rock_channel_atom;

  double  lattice_mag = domain->lattice->xlattice;

  partner0 = NULL;
  partner1 = NULL;
  
  if (atom->nmax > nmax) {
    
    memory->destroy(partner0);
    memory->destroy(partner1);
    nmax = atom->nmax;
    memory->create(partner0,nmax,"fix_injection_update:partner");
    memory->create(partner1,nmax,"fix_injection_update:partner");
    
  }
  
  for (i = 0; i < nall; i++) {
    partner0[i] = FAKE_INT_VALUE;
    partner1[i] = FAKE_INT_VALUE;
  }
  
  for (n = 0; n < nlocal; n++){
    atomi = n;
    if ( atype[n] != CHANNEL_ATOM_TYPE ) continue; // this is not a channel atom
    xi = x0[atomi][0];
    yi = x0[atomi][1];
    zi = x0[atomi][2];
   
    for (m =0; m <nlocal+atom->nghost ; m++){
      atomj =m;
      if ( atype[m] != ROCK_ATOM_TYPE ) continue ; // this is not a rock atom;
      xj = x0[atomj][0];
      yj = x0[atomj][1];
      zj = x0[atomj][2];
      

      if (fabs(fmod(xi,lattice_mag)-0.5*lattice_mag) < 1.e-3*lattice_mag){
	cond0 = ( (fabs(xj - (xi - 0.5*lattice_mag))<1.e-3*lattice_mag) && (fabs(yj - yi) < 1.e-3*lattice_mag) && (fabs(zj - zi)<1.e-3*lattice_mag) );  // first neighbor 
	cond1 = ( (fabs(xj - (xi + 0.5*lattice_mag))<1.e-3*lattice_mag) && (fabs(yj - yi) < 1.e-3*lattice_mag) && (fabs(zj - zi)<1.e-3*lattice_mag) );  
      }
      else if (fabs(fmod(yi,lattice_mag)-0.5*lattice_mag) < 1.e-3*lattice_mag){
	cond0 = ( (fabs(xj - xi)<1.e-3*lattice_mag) && (fabs(yj - (yi - 0.5*lattice_mag)) < 1.e-3*lattice_mag) && (fabs(zj - zi)<1.e-3*lattice_mag) );  
	cond1 = ( (fabs(xj - xi)<1.e-3*lattice_mag) && (fabs(yj - (yi + 0.5*lattice_mag)) < 1.e-3*lattice_mag) && (fabs(zj - zi)<1.e-3*lattice_mag) ); 
      }
      else if (fabs(fmod(zi,lattice_mag)-0.5*lattice_mag) < 1.e-3*lattice_mag){
	cond0 = ( (fabs(xj - xi)<1.e-3*lattice_mag) && (fabs(yj - yi) < 1.e-3*lattice_mag) && (fabs(zj -(zi - 0.5*lattice_mag)) <1.e-3*lattice_mag) );  
	cond1 = ( (fabs(xj - xi)<1.e-3*lattice_mag) && (fabs(yj - yi) < 1.e-3*lattice_mag) && (fabs(zj -(zi + 0.5*lattice_mag)) <1.e-3*lattice_mag) ); 
      }

      
      if (cond0 ==1){partner0[atomi]=tag[atomj];}
      if (cond1 ==1){partner1[atomi]=tag[atomj];}
      
    }
    
  }
  
  


  for (n = 0; n < nlocal; n++){
    atomi = n;

    if (atype[atomi] != CHANNEL_ATOM_TYPE ) continue;
    num_bond_rock_channel_atom[atomi] = 2;
    //    fprintf(screen, "=======bond_rock_channel_atom1 and atom2 %d %d \n",partner0[atomi],partner1[atomi]);
    bond_rock_channel_atom[atomi][0] = partner0[atomi];
    bond_rock_channel_atom[atomi][1] = partner1[atomi];
    
    
  }
  
  //  comm->forward_comm();

}


/*------------------------------------------------------*/
int BondChannelRock::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
 int i,j,m;

  m = 0;
  
  
  for (i = 0; i < n; i++) {
    j = list[i];
    fprintf(screen,"entering pack_comm_fix\n");
    buf[m++] = partner0[j];
    buf[m++] = partner1[j];

  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void BondChannelRock::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){

    partner0[i] = static_cast<int> (buf[m++]);
    partner1[i] = static_cast<int> (buf[m++]);

  }
}
/* ---------------------------------------------------------------------- */

int BondChannelRock::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){

    buf[m++] = partner0[i];
    buf[m++] = partner1[i];
  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void BondChannelRock::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
   
    partner0[j] = static_cast<int> (buf[m++]);
    partner1[j] = static_cast<int> (buf[m++]);
  }
}


double BondChannelRock::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 2*nmax*sizeof(int);
  return bytes;
}






