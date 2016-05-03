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

#include "math.h"
#include "string.h"
#include "compute_channel_pressure.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "bond.h"
#include "memory.h"
#include "error.h"

#include "types.h"

using namespace LAMMPS_NS;

#define DELTA 10000

enum{TIME,ATOMX,ATOMY,ATOMZ,PRESSURE,WIDTH};

/* ---------------------------------------------------------------------- */

ComputeChannelPressure::ComputeChannelPressure(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute bond/local command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute bond/local used when bonds are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  bstyle = new int[nvalues];

  nvalues = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"time") == 0) bstyle[nvalues++] = TIME;
    else if (strcmp(arg[iarg],"atomx") == 0) bstyle[nvalues++] = ATOMX;
    else if (strcmp(arg[iarg],"atomy") == 0) bstyle[nvalues++] = ATOMY;
    else if (strcmp(arg[iarg],"atomz") == 0) bstyle[nvalues++] = ATOMZ;
    else if (strcmp(arg[iarg],"pressure") == 0) bstyle[nvalues++] = PRESSURE;
    else if (strcmp(arg[iarg],"width") == 0) bstyle[nvalues++] = WIDTH;
    else error->all(FLERR,"Invalid keyword in compute bond/local command");
  }

  // set singleflag if need to call bond->single()

  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    //    if (bstyle[i] != DIST) singleflag = 1;

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeChannelPressure::~ComputeChannelPressure()
{
  memory->destroy(vector);
  memory->destroy(array);
  delete [] bstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeChannelPressure::init()
{
  if (force->bond == NULL)
    error->all(FLERR,"No bond style is defined for compute bond/local");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeChannelPressure::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute bond info

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_bonds(1);
}

/* ----------------------------------------------------------------------
   count bonds and compute bond info on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted (type = 0), do not count
   if bond is turned off (type < 0), still count
   if flag is set, compute requested info about bond
   if bond is turned off (type < 0), energy = 0.0
------------------------------------------------------------------------- */

int ComputeChannelPressure::compute_bonds(int flag)
{
  /*
  int i,m,n,atom1,atom2;
  double delx,dely,delz,rsq;
  double *dbuf,*ebuf,*fbuf;
 

  double **x = atom->x;
  int *num_bond = atom->num_bond;
  int *num_bond0 = atom->num_bond0;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *tag = atom->tag;
  int *mask = atom->mask;

  int newton_bond = force->newton_bond;

  Bond *bond = force->bond;
  double eng,fbond;
  double bondx, bondy;
  double **x0 = atom->x0;
  int temp_atom1, temp_atom2;
  double  temp_atom1x,temp_atom1y, temp_atom2x,temp_atom2y;

  //  char buff[50];
  //  int arrsize;
  int bond_per_atom=atom->bond_per_atom;
  m = n = 0;
*/
  int n,m;
  int   *atype  = atom->type;
  //  int *num_bond_rock_channel_atom = atom-> num_bond_rock_channel_atom;
  //  int **bond_rock_channel_atom = atom -> bond_rock_channel_atom;
  double **x0 = atom ->x0;
  double temp_atomx, temp_atomy, temp_atomz;
  double temp_pressure, temp_width;
  double *channel_pressure = atom-> channel_pressure;
  double *channel_width = atom-> channel_width;
  int atom1;
  double *ptr;
  int nlocal = atom->nlocal;
  double time = update->ntimestep *update->dt;
  m = 0 ;
  for (atom1 = 0 ; atom1<nlocal; atom1++){
    if (atype[atom1] != CONNECTED_CHANNEL_ATOM_TYPE) continue;// this is not a wet channel
    //    fprintf(screen, "channel atom %d, x, y ,z %f %f %f p %f width %f \n",atom1,x0[atom1][0],x0[atom1][1],x0[atom1][2], channel_pressure[atom1],channel_width[atom1]);
    temp_atomx = x0[atom1][0];
    temp_atomy = x0[atom1][1];
    temp_atomz = x0[atom1][2];
    temp_pressure = channel_pressure[atom1]*1.e-6; // MPa
    temp_width = channel_width[atom1]*1.e6; //microns

    if (flag){
      
      if (nvalues == 1) ptr = &vector[m];
      else ptr = array[m];
      

      for (n = 0; n < nvalues; n++) {
	//	fprintf(screen,"ptr %f \n", ptr[n]);
	switch (bstyle[n]) {
	case TIME://BONDX:
	  ptr[n] =time;//(double)atom1;//bondx;
	  break;
	case ATOMX://BONDX:
	  ptr[n] =temp_atomx;//(double)atom1;//bondx;
	  break;
	case ATOMY://BONDX:
	  ptr[n] =temp_atomy;//(double)atom1;//bondx;
	  break;
	case ATOMZ://BONDX:
	  ptr[n] =temp_atomz;//(double)atom1;//bondx;
	  break;
	case PRESSURE://BONDX:
	  ptr[n] =temp_pressure;//(double)atom1;//bondx;
	  break;
	case WIDTH://BONDX:
	  ptr[n] =temp_width;//(double)atom1;//bondx;
	  break;
	  
	}
      }
    }
    
    m++;
      
  }
  
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeChannelPressure::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"bond/local:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"bond/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeChannelPressure::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
