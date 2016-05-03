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
#include "check_list.h"
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
#include "domain.h"
#include "irregular.h"

#define TINY  1.e-3 ;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CheckList::CheckList(LAMMPS *lmp) : Pointers(lmp) {}

void CheckList::command(int narg, char **arg)
{
  comm->forward_comm();
  
  if (narg < 1) error->all(FLERR,"Illegal injection_ini command"); // ID group-ID fixchannelini

  check = atof(arg[0]);

  int n;
  int nlocal = atom->nlocal;
  double **x1 = atom->x1;
  double **x = atom->x;
  int *atype = atom->type;
  int *num_bond_channel_atom = atom->num_bond_channel_atom;
  int **bond_channel_atom = atom->bond_channel_atom;
  int temp_atom;
  int max_bond_num, ii;
  double **x0 = atom->x0;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  double temp_x, temp_y, temp_z;
  int *mask = atom->mask;
  int *tag = atom->tag;
  /*
  for (n =0; n<nlocal; n++){
    if (tag[n]==79) fprintf(screen, "=====x,y,z %.2f %.2f %.2f and num_bond %d atype %d, bond_atom1 atom2 atom3 atom4 atom5 %d %d %d %d %d \n",x0[n][0],x0[n][1],x0[n][2],num_bond[n],atype[n],bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3],bond_atom[n][4]);

    if (tag[n]==79) {
      int atom0,atom1, atom2, atom3,atom4;
      atom0 = atom->map(bond_atom[n][0]);
      atom1 = atom->map(bond_atom[n][1]);
      atom2 = atom->map(bond_atom[n][2]);
      atom3 = atom->map(bond_atom[n][3]);
      atom4 = atom->map(bond_atom[n][4]);
    
      fprintf(screen, "=====bond_atom1xyz, 2xyz, 3xyz,4xyz 5xyz(%.0f, %.0f, %.0f), (%.0f, %.0f, %.0f), (%.0f, %.0f, %.0f), (%.0f, %.0f, %.0f), (%.0f, %.0f, %.0f) \n",x0[atom0][0],x0[atom0][1],x0[atom0][2],x0[atom1][0],x0[atom1][1],x0[atom1][2],x0[atom2][0],x0[atom2][1],x0[atom2][2],x0[atom3][0],x0[atom3][1],x0[atom3][2],x0[atom4][0],x0[atom4][1],x0[atom4][2]);
    }
  }
  */

  if (check == 0) {// check rock rock list


    /*
  for (n = 0; n < nlocal; n++){
    if (tag[n] ==681){fprintf(screen, "c==bond_rock1,rock2, rock3, rock4 %d %d %d %d\n",bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3]); }
  }
    */


  
    for (n =0; n<nlocal; n++){
      if (atype[n] != 1) continue;
      x[n][0] = 0.0;
      x[n][1] = 0.0;
      x[n][2] = 0.0;
      max_bond_num = num_bond[n];
      //      if (tag[n] ==79){ fprintf(screen, "==max_bond_num %d, bond_rock1,rock2, rock3, rock4 %d %d %d %d\n",max_bond_num,bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3]); }
      
      for (ii = 0; ii <max_bond_num; ii++){
	
	temp_atom = bond_atom[n][ii];
	temp_atom= atom->map(temp_atom);
	temp_x = x0[temp_atom][0];
	temp_y = x0[temp_atom][1];
	temp_z = x0[temp_atom][2];

	
	//	if (temp_atom <= 0){ temp_x = 0.0; temp_y = 0.0;  temp_z = 0.0;	}
	//	else{
	  //	}

	//	if (n ==660) fprintf(screen, "=====n %d,tag %d, x,y %f %f,neighbor x,y, %f %f,max_bond_num %d  \n",n,tag[n],x0[n][0],x0[n][1],temp_x, temp_y,max_bond_num);
	x[n][0] +=  temp_x/float(max_bond_num);
	x[n][1] +=  temp_y/float(max_bond_num);
	x[n][2] +=  temp_z/float(max_bond_num);
	}
      //      if (tag[n] ==79) fprintf(screen, "n %d, x,y,z %f %f %f \n",tag[n], x[n][0],x[n][1],x[n][2]);
      //      if (tag[n] ==56) fprintf(screen, "n %d, x,y,z %f %f %f \n",tag[n], x[n][0],x[n][1],x[n][2]);
    }
  }

  if (check == 1) {// check channel channel list 
  
    for (n =0; n<nlocal; n++){
      if (atype[n] != 2) continue;
      x[n][0] = 0.0;
      x[n][1] = 0.0;
      x[n][2] = 0.0;
      max_bond_num = num_bond_channel_atom[n];
      //      if (max_bond_num !=12) continue;
      fprintf(screen, "temp_atom %d %d %d %d %d %d %d %d %d %d %d %d max_num_bond %d\n",\
	      bond_channel_atom[n][0],bond_channel_atom[n][1],bond_channel_atom[n][2],bond_channel_atom[n][3],bond_channel_atom[n][4],bond_channel_atom[n][5],\
	      bond_channel_atom[n][6],bond_channel_atom[n][7],bond_channel_atom[n][8],bond_channel_atom[n][9],bond_channel_atom[n][10],bond_channel_atom[n][11],\
	      max_bond_num );

      for (ii = 0; ii <max_bond_num; ii++){

	temp_atom = bond_channel_atom[n][ii];

	if (temp_atom < 0){
	  temp_x = 0.0;
	  temp_y = 0.0;
	  temp_z = 0.0;
	}
	else{
	  temp_atom= atom->map(temp_atom);
	  temp_x = x0[temp_atom][0];
	  temp_y = x0[temp_atom][1];
	  temp_z = x0[temp_atom][2];
	}
	
	x[n][0] +=  temp_x/float(max_bond_num);
	x[n][1] +=  temp_y/float(max_bond_num);
	x[n][2] +=  temp_z/float(max_bond_num);
      }

      fprintf(screen, "x,y,z %f %f %f \n",x[n][0],x[n][1],x[n][2]);
      
    }
  }



  if (check == 2) {// check rock channel list 
    int *num_bond_rock_channel_atom = atom->num_bond_rock_channel_atom;
    int **bond_rock_channel_atom = atom->bond_rock_channel_atom;

    for (n =0; n<nlocal; n++){
      if (atype[n] != 1) continue; // this is not a rock atom
      x[n][0] = 0.0;
      x[n][1] = 0.0;
      x[n][2] = 0.0;
      max_bond_num = num_bond_rock_channel_atom[n];

      for (ii = 0; ii <max_bond_num; ii++){
	
	temp_atom = bond_rock_channel_atom[n][ii];

	if (temp_atom < 0){
	  temp_x = 0.0;
	  temp_y = 0.0;
	  temp_z = 0.0;
	}
	else{
	  temp_atom= atom->map(temp_atom);
	  temp_x = x0[temp_atom][0];
	  temp_y = x0[temp_atom][1];
	  temp_z = x0[temp_atom][2];
	}
	
	x[n][0] +=  temp_x/float(max_bond_num);
	x[n][1] +=  temp_y/float(max_bond_num);
	x[n][2] +=  temp_z/float(max_bond_num);
      }

      fprintf(screen, "x,y,z %f %f %f \n",x[n][0],x[n][1],x[n][2]);
      
    }
  }

  double factor = 1.01;

 if (check == 3) {// check channel rock list 
    int *num_bond_rock_channel_atom = atom->num_bond_rock_channel_atom;
    int **bond_rock_channel_atom = atom->bond_rock_channel_atom;

    for (n =0; n<nlocal; n++){
      if (atype[n] == 1) continue;
      x[n][0] = 0.0;
      x[n][1] = 0.0;
      x[n][2] = 0.0;
      max_bond_num = num_bond_rock_channel_atom[n];
      
      for (ii = 0; ii <max_bond_num; ii++){
	
	temp_atom = bond_rock_channel_atom[n][ii];
	fprintf(screen, "=======max_bond_num %d   rock_atom %d, tag_rock_atom %d \n",max_bond_num, atom->map(temp_atom), temp_atom);
	
	if (temp_atom < 0){
	  temp_x = 0.0;
	  temp_y = 0.0;
	  temp_z = 0.0;
	}
	else{
	  fprintf(screen, "max_bond_num %d   rock_atom %d, tag_rock_atom %d \n",max_bond_num, atom->map(temp_atom), temp_atom);
	  temp_atom= atom->map(temp_atom);
	  temp_x = x0[temp_atom][0];
	  temp_y = x0[temp_atom][1];
	  temp_z = x0[temp_atom][2];
	}
	
	x[n][0] +=  temp_x/float(max_bond_num)*factor;
	x[n][1] +=  temp_y/float(max_bond_num)*factor;
	x[n][2] +=  temp_z/float(max_bond_num)*factor;
      }
      
      fprintf(screen, "x,y,z %f %f %f \n",x[n][0],x[n][1],x[n][2]);
    }
   
 }
 


}

