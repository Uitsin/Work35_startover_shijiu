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
#include "bond_rock_atom.h"
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
#include "irregular.h"
#define TINY  1.e-3 ;
#define FAKE_INT_VALUE  -991 ;

#include "types.h"

//const int ROCK_ATOM_TYPE = 1;
//const int SPRING_BOND_TYPE = 1;

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

BondRockAtom::BondRockAtom(LAMMPS *lmp) : Pointers(lmp) {}

void BondRockAtom::command(int narg, char **arg)
{

  //  if (narg < 3) error->all(FLERR,"Illegal injection_ini command"); // ID group-ID fixchannelini
 

  int nmax =0;
  int nall = atom->nlocal + atom->nghost;

  int n, m,i;
  int *atype = atom->type;
  int *tag = atom->tag;
  int rock_atomi, rock_atomj;
  double xi,yi,zi;
  double xj,yj,zj;
  double **x0 = atom->x0;
  int nlocal = atom->nlocal;
  int Ctype = 0;
  int cond0,cond1,cond2,cond3,cond4,cond5;
  double  lattice_spacing_x=(x0[0][0]-x0[1][0]);
  double  lattice_spacing_y=(x0[0][1]-x0[1][1]);
  double  lattice_mag=sqrt(lattice_spacing_x*lattice_spacing_x+lattice_spacing_y*lattice_spacing_y); //the distance of two neighboring channels

  
  int **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **bond_type =atom->bond_type;
  int ii;
  int max_bond_num;



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


  
  for (n = 0; n < nlocal; n++){
    if ( atype[n] != ROCK_ATOM_TYPE ) continue; // this is not a rock atom
    rock_atomi = n;
    num_bond[rock_atomi] = 6; // total four neighboring rock atoms
    /*
    if ( (fabs(x0[rock_atomi][0]-domain->boxlo[0])<0.1) || (fabs(x0[rock_atomi][0]-domain->boxhi[0]) <0.1)){
      num_bond[rock_atomi]--;
    }
    if ( (fabs(x0[rock_atomi][1]-domain->boxlo[1])<0.1) || (fabs(x0[rock_atomi][1]-domain->boxhi[1]) <0.1)){
      num_bond[rock_atomi]--;
    }
    if ( (fabs(x0[rock_atomi][2]-domain->boxlo[2])<0.1) || (fabs(x0[rock_atomi][2]-domain->boxhi[2]) <0.1)){
      num_bond[rock_atomi]--;
    }
    */
    
    
    for (m =0; m <nlocal+atom->nghost; m++){
      if ( atype[m] != ROCK_ATOM_TYPE ) continue; // this is not a channel atom

      rock_atomj =m;

      xi = x0[rock_atomi][0];
      yi = x0[rock_atomi][1];
      zi = x0[rock_atomi][2];
      
      xj = x0[rock_atomj][0];
      yj = x0[rock_atomj][1];
      zj = x0[rock_atomj][2];
      
      cond0 = ( (fabs(xj - (xi - 1.0) )<1.e-6) && (fabs(yj - yi) < 1.e-6) &&  (fabs(zj - zi) < 1.e-6) );  // first neighbor 
      cond1 = ( (fabs(xj - (xi + 1.0) )<1.e-6) && (fabs(yj - yi) < 1.e-6) &&  (fabs(zj - zi) < 1.e-6) );  // second neighbor
      cond2 = ( (fabs(xj - xi)<1.e-6) && (fabs(yj - (yi - 1.0) ) < 1.e-6) &&  (fabs(zj - zi) < 1.e-6)  );
      cond3 = ( (fabs(xj - xi)<1.e-6) && (fabs(yj - (yi + 1.0) ) < 1.e-6) &&  (fabs(zj - zi) < 1.e-6)  );
      cond4 = ( (fabs(xj - xi)<1.e-6) && (fabs(yj - yi)  < 1.e-6) &&  (fabs(zj - (zi - 1.0)) < 1.e-6)  );
      cond5 = ( (fabs(xj - xi)<1.e-6) && (fabs(yj - yi)  < 1.e-6) &&  (fabs(zj - (zi + 1.0)) < 1.e-6)  );
      
      
      if (cond0 ==1){
	partner0[rock_atomi]=tag[rock_atomj];
	  }
      if (cond1 ==1){
	partner1[rock_atomi]=tag[rock_atomj];
      }
      if (cond2 ==1){
	partner2[rock_atomi]=tag[rock_atomj];
      }
      if (cond3 ==1){
	partner3[rock_atomi]=tag[rock_atomj];
      }
      if (cond4 ==1){
	partner4[rock_atomi]=tag[rock_atomj];
      }
      if (cond5 ==1){
	partner5[rock_atomi]=tag[rock_atomj];
      }

    }

  }


    
  for (n = 0; n < nlocal; n++){
    if ( atype[n] != ROCK_ATOM_TYPE ) continue; // this is not a channel atom
    
    rock_atomi = n;

    int s[6] ={partner0[rock_atomi],partner1[rock_atomi],partner2[rock_atomi],partner3[rock_atomi],partner4[rock_atomi],partner5[rock_atomi]};

    int *p =sort(s,6);

    
    bond_atom[rock_atomi][0] = p[0];
    bond_atom[rock_atomi][1] = p[1];
    bond_atom[rock_atomi][2] = p[2];
    bond_atom[rock_atomi][3] = p[3];
    bond_atom[rock_atomi][4] = p[4];
    bond_atom[rock_atomi][5] = p[5];

    for (ii = 1; ii < 6; ii++){
      if ( bond_atom[rock_atomi][ii-1]>=0 && bond_atom[rock_atomi][ii] < 0){
	num_bond[rock_atomi]= ii;
      }
    }

    max_bond_num = num_bond[rock_atomi];
    for (ii = 0; ii< max_bond_num; ii++){
      bond_type[rock_atomi][ii] =  SPRING_BOND_TYPE ;
    }
  }

  //  comm->forward_comm();

}
/*------------------------------------------------------*/
int BondRockAtom::pack_comm(int n, int *list, double *buf,
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
    buf[m++] = partner4[j];
    buf[m++] = partner5[j];

  }
  return 6;
}

/* ---------------------------------------------------------------------- */

void BondRockAtom::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){

    partner0[i] = static_cast<int> (buf[m++]);
    partner1[i] = static_cast<int> (buf[m++]);
    partner2[i] = static_cast<int> (buf[m++]);
    partner3[i] = static_cast<int> (buf[m++]);
    partner4[i] = static_cast<int> (buf[m++]);
    partner5[i] = static_cast<int> (buf[m++]);

  }
}
/* ---------------------------------------------------------------------- */

int BondRockAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){

    buf[m++] = partner0[i];
    buf[m++] = partner1[i];
    buf[m++] = partner2[i];
    buf[m++] = partner3[i];
    buf[m++] = partner4[i];
    buf[m++] = partner5[i];
  }
  return 6;
}

/* ---------------------------------------------------------------------- */

void BondRockAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
   
    partner0[j] = static_cast<int> (buf[m++]);
    partner1[j] = static_cast<int> (buf[m++]);
    partner2[j] = static_cast<int> (buf[m++]);
    partner3[j] = static_cast<int> (buf[m++]);
    partner4[j] = static_cast<int> (buf[m++]);
    partner5[j] = static_cast<int> (buf[m++]);
  }
}


double BondRockAtom::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 6*nmax*sizeof(int);
  return bytes;
}



int * BondRockAtom::sort(int *s, int n){
  

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
