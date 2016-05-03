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

#include "stdlib.h"
#include "atom_vec_ch_bond.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecChBond::AtomVecChBond(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 1;
  bonds_allow = 1;
  mass_type = 1;

  comm_x_only = 1;
  comm_f_only = 1;
  
  size_forward = 3+3+3+3+3+1;// coords[3]+x0+ f0[3] +channel_pressure+channel_width+channel_width_old+channel_delta_cutoff+tag
  size_reverse = 3; // force
  size_border = 7+3+3+3;
  size_velocity = 3;
  //  size_data_atom = 6;  //1:atom-ID 2:molecule-ID 3:atom-type 4:x 5:y 6:z 
  size_data_atom = 9+3+3;  //1:atom-ID 2:molecule-ID 3:atom-type 4:x 5:y 6:z 7:x0 8:y0 9: z0  10:x1 11:y1 12: z1 
  size_data_vel = 4;
  xcol_data = 4;

  atom->molecule_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecChBond::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  x0 = memory->grow(atom->x0,nmax,3,"atom:x0");
  x1 = memory->grow(atom->x1,nmax,3,"atom:x1");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");
  f0 = memory->grow(atom->f0,nmax,3,"atom:f0");


  molecule = memory->grow(atom->molecule,nmax,"atom:molecule");

  nspecial = memory->grow(atom->nspecial,nmax,3,"atom:nspecial");
  special = memory->grow(atom->special,nmax,atom->maxspecial,"atom:special");

  num_bond = memory->grow(atom->num_bond,nmax,"atom:num_bond");
  num_bond_rock_channel_atom = memory->grow(atom->num_bond_rock_channel_atom,nmax,"atom:num_rock_channel_atom");
  num_bond_channel_atom = memory->grow(atom->num_bond_channel_atom,nmax,"atom:num_channel_atom");

  bond_type = memory->grow(atom->bond_type,nmax,atom->bond_per_atom,"atom:bond_type");

  channel_width = memory->grow(atom->channel_width,nmax,"atom:channel_width");   //channel_width
  channel_width_old = memory->grow(atom->channel_width_old,nmax,"atom:channel_width_old");   //channel_width
  channel_pressure = memory->grow(atom->channel_pressure,nmax,"atom:channel_pressure");   //channel_pressure
  channel_delta_cutoff = memory->grow(atom->channel_delta_cutoff,nmax,"atom:channel_delta_cutoff");   //channel_pressure




  bond_atom = memory->grow(atom->bond_atom,nmax,atom->bond_per_atom,"atom:bond_atom");
  bond_rock_channel_atom = memory->grow(atom->bond_rock_channel_atom,nmax,atom->bond_per_atom,"atom:bond_rock_channel_atom");
  bond_channel_atom = memory->grow(atom->bond_channel_atom,nmax,atom->bond_per_atom,"atom:bond_channel_atom");


  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecChBond::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; 
  x0 =atom->x0;
  x1 =atom->x1;
  v = atom->v; 
  f = atom->f;
  f0 = atom->f0;
  
  molecule = atom->molecule;
  nspecial = atom->nspecial; special = atom->special;
  num_bond = atom->num_bond; 
  num_bond_rock_channel_atom = atom->num_bond_rock_channel_atom;
  num_bond_channel_atom = atom->num_bond_channel_atom;
  bond_type = atom->bond_type;
  channel_width = atom->channel_width;
  channel_width_old = atom->channel_width_old;
  channel_pressure = atom->channel_pressure;
  channel_delta_cutoff = atom->channel_delta_cutoff;
  bond_atom = atom->bond_atom;
  bond_rock_channel_atom = atom->bond_rock_channel_atom;
  bond_channel_atom = atom->bond_channel_atom;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecChBond::copy(int i, int j, int delflag)
{
  int k;
  int l;

  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];

  x0[j][0] = x0[i][0];
  x0[j][1] = x0[i][1];
  x0[j][2] = x0[i][2];

  x1[j][0] = x1[i][0];
  x1[j][1] = x1[i][1];
  x1[j][2] = x1[i][2];

  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];


  f0[j][0] = f0[i][0];
  f0[j][1] = f0[i][1];
  f0[j][2] = f0[i][2];


  molecule[j] = molecule[i];

  num_bond[j] = num_bond[i];
  num_bond_rock_channel_atom[j] = num_bond_rock_channel_atom[i];
  num_bond_channel_atom[j] = num_bond_channel_atom[i];
  channel_width[j]=channel_width[i];
  channel_width_old[j]=channel_width_old[i];
  channel_pressure[j]=channel_pressure[i];
  channel_delta_cutoff[j]=channel_delta_cutoff[i];
 
  for (k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
  }
  for (k = 0; k < num_bond_rock_channel_atom[j]; k++) {
    bond_rock_channel_atom[j][k] = bond_rock_channel_atom[i][k];
  }
  for (k = 0; k < num_bond_channel_atom[j]; k++) {
    bond_channel_atom[j][k] = bond_channel_atom[i][k];
  }
  

  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];
  for (k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::pack_comm(int n, int *list, double *buf,
                           int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];

      buf[m++] = x0[j][0];
      buf[m++] = x0[j][1];
      buf[m++] = x0[j][2];

      buf[m++] = f0[j][0];
      buf[m++] = f0[j][1];
      buf[m++] = f0[j][2];
      

      buf[m++] = channel_width[j];
      buf[m++] = channel_width_old[j];
      buf[m++] = channel_pressure[j];
      buf[m++] = channel_delta_cutoff[j];

      buf[m++] = type[j];



    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      fprintf(screen, "dx dy dz \n");
      buf[m++] = x0[j][0] + dx;
      buf[m++] = x0[j][1] + dy;
      buf[m++] = x0[j][2] + dz;

      buf[m++] = f0[j][0];
      buf[m++] = f0[j][1];
      buf[m++] = f0[j][2];


      buf[m++] = channel_width[j];
      buf[m++] = channel_width_old[j];
      buf[m++] = channel_pressure[j];
      buf[m++] = channel_delta_cutoff[j];

      buf[m++] = type[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::pack_comm_vel(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];

      buf[m++] = x0[j][0];
      buf[m++] = x0[j][1];
      buf[m++] = x0[j][2];

      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];

      buf[m++] = f0[j][0];
      buf[m++] = f0[j][1];
      buf[m++] = f0[j][2];


      buf[m++] = channel_width[j];
      buf[m++] = channel_width_old[j];
      buf[m++] = channel_pressure[j];
      buf[m++] = channel_delta_cutoff[j];
      buf[m++] = type[j];

    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;

        buf[m++] = x0[j][0] + dx;
        buf[m++] = x0[j][1] + dy;
        buf[m++] = x0[j][2] + dz;

        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];

	buf[m++] = f0[j][0];
	buf[m++] = f0[j][1];
	buf[m++] = f0[j][2];


	buf[m++] = channel_width[j];
	buf[m++] = channel_width_old[j];
	buf[m++] = channel_pressure[j];
	buf[m++] = channel_delta_cutoff[j];

	buf[m++] = type[j];

      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;

        buf[m++] = x0[j][0] + dx;
        buf[m++] = x0[j][1] + dy;
        buf[m++] = x0[j][2] + dz;


        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }	
	buf[m++] = f0[j][0];
	buf[m++] = f0[j][1];
	buf[m++] = f0[j][2];


	buf[m++] = channel_width[j];
	buf[m++] = channel_width_old[j];
	buf[m++] = channel_pressure[j];
	buf[m++] = channel_delta_cutoff[j];

	buf[m++] = type[j];
      
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecChBond::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];

    x0[i][0] = buf[m++];
    x0[i][1] = buf[m++];
    x0[i][2] = buf[m++];

    f0[i][0] = buf[m++];
    f0[i][1] = buf[m++];
    f0[i][2] = buf[m++];


    channel_width[i] = buf[m++];
    channel_width_old[i] = buf[m++];
    channel_pressure[i] = buf[m++];
    channel_delta_cutoff[i] = buf[m++];
    type[i] =  static_cast<int> (buf[m++]);


  }
}

/* ---------------------------------------------------------------------- */

void AtomVecChBond::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];

    x0[i][0] = buf[m++];
    x0[i][1] = buf[m++];
    x0[i][2] = buf[m++];

    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];

    f0[i][0] = buf[m++];
    f0[i][1] = buf[m++];
    f0[i][2] = buf[m++];


    channel_width[i] = buf[m++];
    channel_width_old[i] = buf[m++];
    channel_pressure[i] = buf[m++];
    channel_delta_cutoff[i] = buf[m++];
    type[i] =  static_cast<int> (buf[m++]);

  }
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecChBond::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::pack_border(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];

      buf[m++] = x0[j][0];
      buf[m++] = x0[j][1];
      buf[m++] = x0[j][2];

      buf[m++] = x1[j][0];
      buf[m++] = x1[j][1];
      buf[m++] = x1[j][2];

      buf[m++] = f0[j][0];
      buf[m++] = f0[j][1];
      buf[m++] = f0[j][2];


      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;

      buf[m++] = x0[j][0] + dx;
      buf[m++] = x0[j][1] + dy;
      buf[m++] = x0[j][2] + dz;

      buf[m++] = x1[j][0] + dx;
      buf[m++] = x1[j][1] + dy;
      buf[m++] = x1[j][2] + dz;

      buf[m++] = f0[j][0];
      buf[m++] = f0[j][1];
      buf[m++] = f0[j][2];

      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::pack_border_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];

      buf[m++] = x0[j][0];
      buf[m++] = x0[j][1];
      buf[m++] = x0[j][2];

      buf[m++] = x1[j][0];
      buf[m++] = x1[j][1];
      buf[m++] = x1[j][2];

      buf[m++] = f0[j][0];
      buf[m++] = f0[j][1];
      buf[m++] = f0[j][2];

      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;

        buf[m++] = x0[j][0] + dx;
        buf[m++] = x0[j][1] + dy;
        buf[m++] = x0[j][2] + dz;

        buf[m++] = x1[j][0] + dx;
        buf[m++] = x1[j][1] + dy;
        buf[m++] = x1[j][2] + dz;

	buf[m++] = f0[j][0];
	buf[m++] = f0[j][1];
	buf[m++] = f0[j][2];

        buf[m++] = tag[j];
        buf[m++] = type[j];
        buf[m++] = mask[j];
        buf[m++] = molecule[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;

        buf[m++] = x0[j][0] + dx;
        buf[m++] = x0[j][1] + dy;
        buf[m++] = x0[j][2] + dz;

        buf[m++] = x1[j][0] + dx;
        buf[m++] = x1[j][1] + dy;
        buf[m++] = x1[j][2] + dz;

	buf[m++] = f0[j][0];
	buf[m++] = f0[j][1];
	buf[m++] = f0[j][2];

        buf[m++] = tag[j];
        buf[m++] = type[j];
        buf[m++] = mask[j];
        buf[m++] = molecule[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = molecule[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecChBond::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];

    x0[i][0] = buf[m++];
    x0[i][1] = buf[m++];
    x0[i][2] = buf[m++];

    x1[i][0] = buf[m++];
    x1[i][1] = buf[m++];
    x1[i][2] = buf[m++];

    f0[i][0] = buf[m++];
    f0[i][1] = buf[m++];
    f0[i][2] = buf[m++];

    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    molecule[i] = static_cast<int> (buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecChBond::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];

    x0[i][0] = buf[m++];
    x0[i][1] = buf[m++];
    x0[i][2] = buf[m++];

    x1[i][0] = buf[m++];
    x1[i][1] = buf[m++];
    x1[i][2] = buf[m++];

    f0[i][0] = buf[m++];
    f0[i][1] = buf[m++];
    f0[i][2] = buf[m++];

    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    molecule[i] = static_cast<int> (buf[m++]);
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    molecule[i] = static_cast<int> (buf[m++]);
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecChBond::pack_exchange(int i, double *buf)
{
  int k;
  int l;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];

  buf[m++] = x0[i][0];
  buf[m++] = x0[i][1];
  buf[m++] = x0[i][2];

  buf[m++] = x1[i][0];
  buf[m++] = x1[i][1];
  buf[m++] = x1[i][2];

  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = f[i][0];
  buf[m++] = f[i][1];
  buf[m++] = f[i][2];


  buf[m++] = f0[i][0];
  buf[m++] = f0[i][1];
  buf[m++] = f0[i][2];

  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  *((tagint *) &buf[m++]) = image[i];

  buf[m++] = molecule[i];

  buf[m++] = num_bond[i];
  buf[m++] = num_bond_rock_channel_atom[i];
  buf[m++] = num_bond_channel_atom[i];

  buf[m++] = channel_width[i];
  buf[m++] = channel_width_old[i];
  buf[m++] = channel_pressure[i];
  buf[m++] = channel_delta_cutoff[i];

  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = bond_type[i][k];
    buf[m++] = bond_atom[i][k];
  }
  for (k = 0; k < num_bond_rock_channel_atom[i]; k++) {
    buf[m++] = bond_rock_channel_atom[i][k];
  }
  for (k = 0; k < num_bond_channel_atom[i]; k++) {
    buf[m++] = bond_channel_atom[i][k];
  }
   


  buf[m++] = nspecial[i][0];
  buf[m++] = nspecial[i][1];
  buf[m++] = nspecial[i][2];
  for (k = 0; k < nspecial[i][2]; k++) buf[m++] = special[i][k];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecChBond::unpack_exchange(double *buf)
{
  int k;
  int l;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];

  x0[nlocal][0] = buf[m++];
  x0[nlocal][1] = buf[m++];
  x0[nlocal][2] = buf[m++];

  x1[nlocal][0] = buf[m++];
  x1[nlocal][1] = buf[m++];
  x1[nlocal][2] = buf[m++];

  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  f[nlocal][0] = buf[m++];
  f[nlocal][1] = buf[m++];
  f[nlocal][2] = buf[m++];


  f0[nlocal][0] = buf[m++];
  f0[nlocal][1] = buf[m++];
  f0[nlocal][2] = buf[m++];

  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = *((tagint *) &buf[m++]);

  molecule[nlocal] = static_cast<int> (buf[m++]);

  num_bond[nlocal] = static_cast<int> (buf[m++]);
  num_bond_rock_channel_atom[nlocal] = static_cast<int> (buf[m++]);
  num_bond_channel_atom[nlocal] = static_cast<int> (buf[m++]);

  channel_width[nlocal] = buf[m++];
  channel_width_old[nlocal] = buf[m++];
  channel_pressure[nlocal] = buf[m++];
  channel_delta_cutoff[nlocal] = buf[m++];

  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = static_cast<int> (buf[m++]);
    bond_atom[nlocal][k] = static_cast<int> (buf[m++]);  
  }
  for (k = 0; k < num_bond_rock_channel_atom[nlocal]; k++) {
    bond_rock_channel_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }
  for (k = 0; k < num_bond_channel_atom[nlocal]; k++) {
    bond_channel_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }


  nspecial[nlocal][0] = static_cast<int> (buf[m++]);
  nspecial[nlocal][1] = static_cast<int> (buf[m++]);
  nspecial[nlocal][2] = static_cast<int> (buf[m++]);
  for (k = 0; k < nspecial[nlocal][2]; k++)
    special[nlocal][k] = static_cast<int> (buf[m++]);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecChBond::size_restart()
{
  int i;
  int k;
  int nlocal = atom->nlocal;
  int n = 0;
  for (i = 0; i < nlocal; i++){
    //    n += 16 + 2*num_bond[i]; //13+3
    //    n += 17 + 3*num_bond0[i]; //13+3+1//13+x0,y0,z0+x1[3]+ f0[3]+ num_bond0 + channel_index +num_bond_rock_channel_atom +num_bond_channel_atom+channel_index[i],channel_neighbor_atom2,channel_neighbor_atom3,channel_neighbor_atom4, channel_width, channel_pressure
    n += 34+ 6*num_bond[i];//6 lists, which include bond_type, channel_type,channel_width, channel_pressure,bond_atom,num_channel_neighbor, bond_rock_channel_atom, bond_channel_atom

  }
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecChBond::pack_restart(int i, double *buf)
{
  int k;
  int l;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];

  buf[m++] = x0[i][0];
  buf[m++] = x0[i][1];
  buf[m++] = x0[i][2];

  buf[m++] = x1[i][0];
  buf[m++] = x1[i][1];
  buf[m++] = x1[i][2];

  buf[m++] = f0[i][0];
  buf[m++] = f0[i][1];
  buf[m++] = f0[i][2];

  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  *((tagint *) &buf[m++]) = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = molecule[i];

  buf[m++] = num_bond[i];
  buf[m++] = num_bond_rock_channel_atom[i];
  buf[m++] = num_bond_channel_atom[i];
  
  buf[m++] = channel_width[i];
  buf[m++] = channel_width_old[i];
  buf[m++] = channel_pressure[i];
  buf[m++] = channel_delta_cutoff[i];

  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = MAX(bond_type[i][k],-bond_type[i][k]);
    buf[m++] = bond_atom[i][k];   
  }
  for (k = 0; k < num_bond_rock_channel_atom[i]; k++) {
    buf[m++] = bond_rock_channel_atom[i][k];
  }
  for (k = 0; k < num_bond_channel_atom[i]; k++) {
    buf[m++] = bond_channel_atom[i][k];
  }
   


  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecChBond::unpack_restart(double *buf)
{
  int k;
  int l;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];

  x0[nlocal][0] = buf[m++];
  x0[nlocal][1] = buf[m++];
  x0[nlocal][2] = buf[m++];

  x1[nlocal][0] = buf[m++];
  x1[nlocal][1] = buf[m++];
  x1[nlocal][2] = buf[m++];

  f0[nlocal][0] = buf[m++];
  f0[nlocal][1] = buf[m++];
  f0[nlocal][2] = buf[m++];

  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = *((tagint *) &buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  molecule[nlocal] = static_cast<int> (buf[m++]);

  num_bond[nlocal] = static_cast<int> (buf[m++]);
  num_bond_rock_channel_atom[nlocal] = static_cast<int> (buf[m++]);
  num_bond_channel_atom[nlocal] = static_cast<int> (buf[m++]);

  channel_width[nlocal] = buf[m++];
  channel_width_old[nlocal] = buf[m++];
  channel_pressure[nlocal] = buf[m++];
  channel_delta_cutoff[nlocal] = buf[m++];

  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = static_cast<int> (buf[m++]);
    bond_atom[nlocal][k] = static_cast<int> (buf[m++]);    
  }
  for (k = 0; k < num_bond_rock_channel_atom[nlocal]; k++) {
    bond_rock_channel_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }
  for (k = 0; k < num_bond_channel_atom[nlocal]; k++) {
    bond_channel_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }


  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecChBond::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  int k;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  x0[nlocal][0] = coord[0];
  x0[nlocal][1] = coord[1];
  x0[nlocal][2] = coord[2];

  x1[nlocal][0] = coord[0];
  x1[nlocal][1] = coord[1];
  x1[nlocal][2] = coord[2];

  f0[nlocal][0] = 0.0;
  f0[nlocal][1] = 0.0;
  f0[nlocal][2] = 0.0;

  mask[nlocal] = 1;
  image[nlocal] = ((tagint) IMGMAX << IMG2BITS) |
    ((tagint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  molecule[nlocal] = 0;
  num_bond[nlocal] = 0;
  num_bond_rock_channel_atom[nlocal] = 0;
  num_bond_channel_atom[nlocal] = 0;

  channel_width[nlocal] =0;
  channel_width_old[nlocal] =0;
  channel_pressure[nlocal] =0;
  channel_delta_cutoff[nlocal] =0;

  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecChBond::data_atom(double *coord, tagint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  int k;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");

  molecule[nlocal] = atoi(values[1]);

  type[nlocal] = atoi(values[2]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  x0[nlocal][0] = coord[0];
  x0[nlocal][1] = coord[1];
  x0[nlocal][2] = coord[2];

  x1[nlocal][0] = coord[0];
  x1[nlocal][1] = coord[1];
  x1[nlocal][2] = coord[2];

  f0[nlocal][0] = 0.0;
  f0[nlocal][1] = 0.0;
  f0[nlocal][2] = 0.0;

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  num_bond[nlocal] = 0;
  num_bond_rock_channel_atom[nlocal] = 0;
  num_bond_channel_atom[nlocal] = 0;

  channel_width[nlocal] =0;
  channel_width_old[nlocal] =0;
  channel_pressure[nlocal] =0;
  channel_delta_cutoff[nlocal] =0;


  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecChBond::data_atom_hybrid(int nlocal, char **values)
{
  molecule[nlocal] = atoi(values[0]);

  num_bond[nlocal] = 0;


  return 1;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecChBond::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = tag[i];
    buf[i][1] = molecule[i];
    buf[i][2] = type[i];
    buf[i][3] = x[i][0];
    buf[i][4] = x[i][1];
    buf[i][5] = x[i][2];

    buf[i][6] = x0[i][0];
    buf[i][7] = x0[i][1];
    buf[i][8] = x0[i][2];

    buf[i][9] = x1[i][0];
    buf[i][10] = x1[i][1];
    buf[i][11] = x1[i][2];

    buf[i][12] = (image[i] & IMGMASK) - IMGMAX;
    buf[i][13] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    buf[i][14] = (image[i] >> IMG2BITS) - IMGMAX;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecChBond::pack_data_hybrid(int i, double *buf)
{
  buf[0] = molecule[i];
  return 1;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecChBond::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d %d  %g %g %g %g %g %g %d %d %d\n",
            (int) buf[i][0],(int) buf[i][1],(int) buf[i][2],
            buf[i][3],buf[i][4],buf[i][5],
            buf[i][6],buf[i][7],buf[i][8],
            buf[i][9],buf[i][10],buf[i][11],
            (int) buf[i][12],(int) buf[i][13],(int) buf[i][14]);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecChBond::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %d",(int) buf[0]);
  return 1;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecChBond::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);

  if (atom->memcheck("x0")) bytes += memory->usage(x0,nmax,3); // adding original position of each atom
  if (atom->memcheck("x1")) bytes += memory->usage(x1,nmax,3); // adding original position of each atom

  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);
  if (atom->memcheck("f0")) bytes += memory->usage(f0,nmax,3);

  if (atom->memcheck("molecule")) bytes += memory->usage(molecule,nmax);
  if (atom->memcheck("nspecial")) bytes += memory->usage(nspecial,nmax,3);
  if (atom->memcheck("special"))
    bytes += memory->usage(special,nmax,atom->maxspecial);

  if (atom->memcheck("num_bond")) bytes += memory->usage(num_bond,nmax);
  if (atom->memcheck("num_bond_rock_channel_atom")) bytes += memory->usage(num_bond_rock_channel_atom,nmax);
  if (atom->memcheck("num_bond_channel_atom")) bytes += memory->usage(num_bond_channel_atom,nmax);

  if (atom->memcheck("bond_type"))
    bytes += memory->usage(bond_type,nmax,atom->bond_per_atom);

  if (atom->memcheck("channel_width")) bytes += memory->usage(channel_width,nmax);
  if (atom->memcheck("channel_width_old")) bytes += memory->usage(channel_width_old,nmax);
  if (atom->memcheck("channel_pressure")) bytes += memory->usage(channel_pressure,nmax);
  if (atom->memcheck("channel_delta_cutoff")) bytes += memory->usage(channel_delta_cutoff,nmax);

  if (atom->memcheck("bond_atom"))bytes += memory->usage(bond_atom,nmax,atom->bond_per_atom);
  if (atom->memcheck("bond_rock_channel_atom"))bytes += memory->usage(bond_rock_channel_atom,nmax,atom->bond_per_atom);
  if (atom->memcheck("bond_channel_atom"))bytes += memory->usage(bond_channel_atom,nmax,atom->bond_per_atom);

  return bytes;
}
