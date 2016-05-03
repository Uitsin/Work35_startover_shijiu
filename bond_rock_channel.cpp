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
#include "bond_rock_channel.h"
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

BondRockChannel::BondRockChannel(LAMMPS *lmp) : Pointers(lmp) {}

void BondRockChannel::command(int narg, char **arg)
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
  int ii;

  double  lattice_mag = domain->lattice->xlattice;

  partner0 = NULL;
  partner1 = NULL;
  partner2 = NULL;
  partner3 = NULL;
  partner4 = NULL;
  partner5 = NULL;
  
  if (atom->nmax > nmax) {
    
    memory->destroy(partner0);
    memory->destroy(partner1);
    memory->destroy(partner2);
    memory->destroy(partner3);
    memory->destroy(partner4);
    memory->destroy(partner5);
    
    nmax = atom->nmax;
    
    memory->create(partner0,nmax,"fix_injection_update:partner");
    memory->create(partner1,nmax,"fix_injection_update:partner");
    memory->create(partner2,nmax,"fix_injection_update:partner");
    memory->create(partner3,nmax,"fix_injection_update:partner");
    memory->create(partner4,nmax,"fix_injection_update:partner");
    memory->create(partner5,nmax,"fix_injection_update:partner");
    
  }
  
  
  for (i = 0; i < nall; i++) {
    partner0[i] = FAKE_INT_VALUE;
    partner1[i] = FAKE_INT_VALUE;
    partner2[i] = FAKE_INT_VALUE;
    partner3[i] = FAKE_INT_VALUE;
    partner4[i] = FAKE_INT_VALUE;
    partner5[i] = FAKE_INT_VALUE;
    
  }
  
  int **bond_atom = atom -> bond_atom;
  //  for (n = 0; n < nlocal; n++){ if (tag[n] ==681){fprintf(screen, "a==bond_rock1,rock2, rock3, rock4 %d %d %d %d\n",bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3]);} }


  
  for (n = 0; n < nlocal; n++){
    atomi = n;
    
    
    xi = x0[atomi][0];
    yi = x0[atomi][1];
    zi = x0[atomi][2];
    
   
    for (m =0; m <nlocal+atom->nghost ; m++){
      atomj =m;
      xj = x0[atomj][0];
      yj = x0[atomj][1];
      zj = x0[atomj][2];
      
      if ( atype[atomi]== ROCK_ATOM_TYPE && atype[atomj] == CHANNEL_ATOM_TYPE ){ 

       	cond0 = ( (fabs(xj - (xi - 0.5*lattice_mag))<1.e-3) && (fabs(yj - yi) < 1.e-3) && (fabs(zj -zi) <1.e-3) );  // first neighbor 
       	cond1 = ( (fabs(xj - (xi + 0.5*lattice_mag))<1.e-3) && (fabs(yj - yi) < 1.e-3) && (fabs(zj -zi) <1.e-3) );  
       	cond2 = ( (fabs(xj - xi)<1.e-3) && (fabs(yj - (yi -0.5*lattice_mag) ) < 1.e-3) && (fabs(zj -zi) <1.e-3) );  
       	cond3 = ( (fabs(xj - xi)<1.e-3) && (fabs(yj - (yi +0.5*lattice_mag) ) < 1.e-3) && (fabs(zj -zi) <1.e-3) ); 
       	cond4 = ( (fabs(xj - xi)<1.e-3) && (fabs(yj - yi) < 1.e-3) && (fabs(zj -(zi-0.5*lattice_mag) ) <1.e-3) );  
       	cond5 = ( (fabs(xj - xi)<1.e-3) && (fabs(yj - yi) < 1.e-3) && (fabs(zj -(zi+0.5*lattice_mag) ) <1.e-3) ); 



	if (cond0 ==1){partner0[atomi]=tag[atomj];}
	if (cond1 ==1){partner1[atomi]=tag[atomj];}
	if (cond2 ==1){partner2[atomi]=tag[atomj];}
	if (cond3 ==1){partner3[atomi]=tag[atomj];}
	if (cond4 ==1){partner4[atomi]=tag[atomj];}
	if (cond5 ==1){partner5[atomi]=tag[atomj];}

      }
      
    }
    
  }




  for (n = 0; n < nlocal; n++){
    atomi = n;


    if (atype[atomi] == ROCK_ATOM_TYPE) { // this is a rock atom
      num_bond_rock_channel_atom[atomi] = 6;
      /*
      if ( fabs(x0[atomi][0]-domain->boxlo[0])<0.1 ){
	num_bond_rock_channel_atom[atomi]--;
      }

      if (fabs(x0[atomi][0]-domain->boxhi[0]) <0.1 ) {
	num_bond_rock_channel_atom[atomi]--;
      }
      if ( fabs(x0[atomi][1]-domain->boxlo[1])<0.1 ){
	num_bond_rock_channel_atom[atomi]--;
      }

      if (fabs(x0[atomi][1]-domain->boxhi[1]) <0.1 ) {
	num_bond_rock_channel_atom[atomi]--;
      }

      if ( fabs(x0[atomi][2]-domain->boxlo[2])<0.1 ){
	num_bond_rock_channel_atom[atomi]--;
      }

      if (fabs(x0[atomi][2]-domain->boxhi[2]) <0.1 ) {
	num_bond_rock_channel_atom[atomi]--;
      }
      */

      int s[6] ={partner0[atomi],partner1[atomi],partner2[atomi],partner3[atomi],partner4[atomi],partner5[atomi]};
      int *p =sort(s,6);

  
      bond_rock_channel_atom[atomi][0] = p[0];
      bond_rock_channel_atom[atomi][1] = p[1];
      bond_rock_channel_atom[atomi][2] = p[2]; 
      bond_rock_channel_atom[atomi][3] = p[3];
      bond_rock_channel_atom[atomi][4] = p[4];
      bond_rock_channel_atom[atomi][5] = p[5];



      
      for (ii = 1; ii < 6; ii++){
	if ( bond_rock_channel_atom[atomi][ii-1]>=0 && bond_rock_channel_atom[atomi][ii] < 0){
	  num_bond_rock_channel_atom[atomi] = ii ;
	}
      }
      
      
    }
    
  }
  
  //  comm->forward_comm();


}


/*
int BondRockChannel::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
 int i,j,m;

  m = 0;
  
  
  for (i = 0; i < n; i++) {
    j = list[i];

    buf[m++] = partner0[j];
    buf[m++] = partner1[j];
    buf[m++] = partner2[j];
    buf[m++] = partner3[j];

  }
  return 4;
}


void BondRockChannel::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){

    partner0[i] = static_cast<int> (buf[m++]);
    partner1[i] = static_cast<int> (buf[m++]);
    partner2[i] = static_cast<int> (buf[m++]);
    partner3[i] = static_cast<int> (buf[m++]);

  }
}

int BondRockChannel::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){

    buf[m++] = partner0[i];
    buf[m++] = partner1[i];
    buf[m++] = partner2[i];
    buf[m++] = partner3[i];
  }
  return 4;
}



void BondRockChannel::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
   
    partner0[j] = static_cast<int> (buf[m++]);
    partner1[j] = static_cast<int> (buf[m++]);
    partner2[j] = static_cast<int> (buf[m++]);
    partner3[j] = static_cast<int> (buf[m++]);
  }
}


double BondRockChannel::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 4*nmax*sizeof(int);
  return bytes;
}

*/


int * BondRockChannel::sort(int *s, int n){
  

  int ii, jj,temp;

  for(int ii=0; ii<n; ii++){
    for(int jj=0; jj<n-1; jj++){
      if(s[jj]<s[jj+1]){
	temp = s[jj+1];
	s[jj+1] = s[jj];
	s[jj] = temp;
      }
    }
  }
  return s;
}


