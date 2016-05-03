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
#include <iostream>     // std::cout
#include <limits>       // std::numeric_limits

#include "string.h"
#include "math.h"
#include "fix_lubrication_update_BC.h"
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
#define  FAKE_INT_VALUE -991;
using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixLubricationUpdateBC::FixLubricationUpdateBC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal fix lubrication/update command"); // ID group-ID fixchannelini

  BC_xlo = atof(arg[3]);
  BC_xhi = atof(arg[4]);

  BC_ylo = atof(arg[5]);
  BC_yhi = atof(arg[6]);

  BC_zlo = atof(arg[7]);
  BC_zhi = atof(arg[8]);

  vij_max = atof(arg[9]);

  partner = NULL;
  nmax =0;
}



FixLubricationUpdateBC::~FixLubricationUpdateBC()
{
  memory->destroy(partner);
}


/*============
  setmask
  ============*/
int FixLubricationUpdateBC::setmask()
{
  int mask = 0;
  //   mask |= PRE_FORCE;
   mask |= FINAL_INTEGRATE;
  return mask;
}


/*==================================
  called before force routine
  ==================================*/

//void FixLubricationUpdateBC::pre_force(int vflag)
void FixLubricationUpdateBC::final_integrate()

{
  
  if (atom->nmax > nmax) {
    memory->destroy(partner);
    nmax = atom->nmax;
    memory->create(partner,nmax,"fix_lubrication_update:partner");
  }

  //comm->borders();
  comm->forward_comm();
  new_channel();
  comm->forward_comm();

  bond_break();
  neighbor->build_topology();
  comm->forward_comm();


  comm->forward_comm();
  cal_channel_pressure();
  comm->forward_comm();

  copy_channel_width();
  comm->forward_comm();

  lubrication();
  comm->forward_comm();

  neighbor->build_topology();
  channel_update();
  comm->forward_comm();



}


/*==================================
  copy channel_width to channel_width_old
  ==================================*/

void FixLubricationUpdateBC::copy_channel_width()
{
  int nlocal = atom->nlocal;
  int *atype = atom->type;
  double *channel_width = atom->channel_width;
  double *channel_width_old = atom->channel_width_old;
  int n;
  for (n=0; n<nlocal;n++){
    if (atype[n] != CONNECTED_CHANNEL_ATOM_TYPE) continue;
    channel_width_old[n]= channel_width[n];
  }
}

/*==================================
  calculate channel pressure
  ==================================*/
void FixLubricationUpdateBC::cal_channel_pressure()
{
  
  int n;
  int nlocal = atom->nlocal;
  int *atype = atom->type;
  int *num_bond_rock_channel_atom = atom->num_bond_rock_channel_atom;
  int **bond_rock_channel_atom = atom->bond_rock_channel_atom;
  int rock_atom1, rock_atom2,channel_atom;
  double **x0 = atom->x0;
  double **x = atom->x;
  double **f0 = atom->f0;
  double delta_x,delta_y,delta_z;
  double pressure;
  double *channel_pressure = atom->channel_pressure;
  double *channel_width = atom->channel_width;
  double delta_mag;
  int *tag = atom->tag;
  double check_force;
  for (n=0; n<nlocal; n++){
    if (atype[n] != CONNECTED_CHANNEL_ATOM_TYPE) continue;
    channel_atom = n;
    rock_atom1 = atom->map(bond_rock_channel_atom[channel_atom][0]);
    rock_atom2 = atom->map(bond_rock_channel_atom[channel_atom][1]);
    
    delta_x = (x0[rock_atom1][0] - x0[rock_atom2][0]);
    delta_y = (x0[rock_atom1][1] - x0[rock_atom2][1]);
    delta_z = (x0[rock_atom1][2] - x0[rock_atom2][2]);
    
    delta_mag = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
    check_force = (f0[rock_atom1][0]*delta_x +  f0[rock_atom1][1]*delta_y + f0[rock_atom1][2]*delta_z )/delta_mag;

    if ( (update->ntimestep > 10)  && check_force ==0 ) fprintf(screen, "Wrong!! force is zero !!!%f at timestep %d\n", check_force, update->ntimestep);

    pressure =-( (f0[rock_atom1][0]*delta_x+f0[rock_atom1][1]*delta_y + f0[rock_atom1][2]*delta_z) -(f0[rock_atom2][0]*delta_x+f0[rock_atom2][1]*delta_y + f0[rock_atom2][2]*delta_z) )/2.0/delta_mag;

    if (pressure < 0) {
      //      fprintf(screen, "at t= %.3f, this channel has negative pressure! %.2f (MPa) at x y z %f %f %f \n", update->ntimestep*update->dt, pressure*1.e-6, x0[channel_atom][0],x0[channel_atom][1],x0[channel_atom][2]);
      pressure = 0; // this channel is partially filled with water
    }

    channel_pressure[channel_atom] = pressure;

      if ((update->ntimestep %50000)==0) fprintf(screen, "channel pressure is %.2f MPa at x y z %f %f %f  \n",channel_pressure[channel_atom]*1.e-6,x0[channel_atom][0],x0[channel_atom][1],x0[channel_atom][2]);
    
  }

}



/*========================================
  channel flow based on lubrication approximation
  ========================================*/

void FixLubricationUpdateBC::lubrication(){

  int n;
  int nlocal = atom->nlocal;
  int *atype = atom->type;
  double **x0 = atom ->x0;
  double **x = atom ->x;
  double channel_wi, channel_wj;
  double channel_pi, channel_pj;
  double local_w;
  double on_twelve_mu_L;
  double **v = atom->v;
  double L = 1.0;
  double vil; 
  double mu = 1.e-3; //Pa.s
  int channel_atomi, channel_atomj;
  int *num_bond_channel_atom = atom ->num_bond_channel_atom;
  int **bond_channel_atom = atom->bond_channel_atom;
  int **bond_rock_channel_atom = atom->bond_rock_channel_atom;

  int ii;
  int *tag = atom->tag;
  double *channel_width = atom->channel_width;
  double *channel_width_old = atom->channel_width_old;
  double *channel_pressure = atom->channel_pressure;
  double vij;

  int nstep = update->ntimestep;
  double ntime = update->dt;
  double dw;
  on_twelve_mu_L =1.0/(12.0*mu*L);// need to be rescaled if lattice spacing L is not 1 meter.

  for (n=0; n<nlocal;n++){
    if (atype[n] != CONNECTED_CHANNEL_ATOM_TYPE ) continue;
    channel_atomi = n;
    for (ii = 0; ii <num_bond_channel_atom[channel_atomi]; ii++){
      //      if (tag[channel_atomi] > tag[channel_atomj]) continue; // double counting 
      channel_atomj = atom->map(bond_channel_atom[channel_atomi][ii] );
      if (atype[channel_atomj] != CONNECTED_CHANNEL_ATOM_TYPE ) continue;
      channel_wi = channel_width_old[channel_atomi];
      channel_pi = channel_pressure[channel_atomi];
      
      channel_wj = channel_width_old[channel_atomj];
      channel_pj = channel_pressure[channel_atomj];

      local_w =0; 

      if (channel_pi > channel_pj){local_w = channel_wi;}
      else {local_w = channel_wj;}
      vij = (channel_pj-channel_pi)*local_w*local_w*local_w*on_twelve_mu_L;
      if ((channel_pi < 0) || (channel_pj < 0)){fprintf(screen, "Negative pressure!!! pressure calculation is wrong!!!!!!. Check again!! \n");}
      if (fabs(vij) > vij_max){
	if (vij > 0){vij = vij_max;}
	else {vij = -vij_max;}  
      }

      dw = vij * update->dt;
      channel_width[channel_atomi]+=dw;
    }
  }
}


/*========================================
  channel atom update based on new channel width
  ========================================*/
void FixLubricationUpdateBC::channel_update(){

  int n;
  int nbondlist = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;
  int rock_atom1, rock_atom2,channel_atom;
  int btype;
  int *atype = atom->type;
  double  *channel_width = atom->channel_width;
  double  *channel_width_old = atom->channel_width_old;
  double channel_w,channel_w0;
  double delta_x, delta_y, delta_z, delta_mag;
  double **x0 = atom->x0;
  double **x = atom->x;
  int newton_bond = force->newton_bond;
  int nlocal = atom->nlocal;

  //  neighbor->build_topology();

  for (n = 0; n < nbondlist; n++) {
    rock_atom1 = bondlist[n][0];
    rock_atom2 = bondlist[n][1];
    btype = bondlist[n][2];
    if (btype != WET_CHANNEL_BOND_TYPE) continue;// this is not a wet channel

    channel_atom = find_channel_atom(rock_atom1, rock_atom2);

    if (atype[channel_atom] != CONNECTED_CHANNEL_ATOM_TYPE ) {fprintf(screen, "!!!!this bond should be wet channel but it is %d. Check again \n",atype[channel_atom]); continue;}

    channel_w = channel_width[channel_atom];
      
    delta_x = x0[rock_atom1][0] - x0[rock_atom2][0];
    delta_y = x0[rock_atom1][1] - x0[rock_atom2][1];
    delta_z = x0[rock_atom1][2] - x0[rock_atom2][2];
    delta_mag = sqrt (delta_x *delta_x + delta_y *delta_y +delta_z *delta_z);
    

    if (newton_bond || rock_atom1 <nlocal){
      if (fabs(delta_x)  == delta_mag) {
	x[rock_atom1][0] = x[channel_atom][0] +0.5*(1+channel_w) * delta_x/delta_mag;
      }
      if (fabs(delta_y)  == delta_mag) {
	x[rock_atom1][1] = x[channel_atom][1] +0.5*(1+channel_w) * delta_y/delta_mag;
      }
      if (fabs(delta_z)  == delta_mag) {
	x[rock_atom1][2] = x[channel_atom][2] +0.5*(1+channel_w) * delta_z/delta_mag;
      }

    }
    if (newton_bond || rock_atom2 <nlocal){
      if (fabs(delta_x)  == delta_mag) {
	x[rock_atom2][0] = x[channel_atom][0] -0.5*(1+channel_w) * delta_x/delta_mag;
      }
      if (fabs(delta_y)  == delta_mag) {
	x[rock_atom2][1] = x[channel_atom][1] -0.5*(1+channel_w) * delta_y/delta_mag;
      }
      if (fabs(delta_z)  == delta_mag) {
	x[rock_atom2][2] = x[channel_atom][2] -0.5*(1+channel_w) * delta_z/delta_mag;
      }
    }
  }
}




/*====================================================
  find channel_atom between rock_atom1, and rock_atom2
  ====================================================*/

int FixLubricationUpdateBC::find_channel_atom(int rock_atom1, int rock_atom2){
  
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



/*==================================
if a ISOLATED_CHANNEL_ATOM_TYPE atom is now a CONNECTED_CHANNEL_ATOM_TYPE atom, it becomes CONNECTED_CHANNEL_ATOM_TYPE.  
  ==================================*/

void FixLubricationUpdateBC::new_channel()
{ 
  int *atype = atom->type;
  int n;
  int nlocal = atom->nlocal;
  int channel_atom;
  int rock_atom1, rock_atom2;
  double delta_x, delta_y, delta_z,dist;
  double **x = atom->x;
  int **bond_rock_channel_atom = atom->bond_rock_channel_atom;
  double *channel_width = atom->channel_width;
  int max_bond_num; 
  int *num_bond_channel_atom = atom->num_bond_channel_atom;
  int **bond_channel_atom = atom->bond_channel_atom;
  int ii;
  int channel_atomj;
  int connection_check[6];
  
  int check;
  int m;
  int cond0,cond1,cond2, cond3, cond4,cond5;
  int temp_channel_atom;
  double **x0 = atom->x0;
  double *channel_delta_cutoff = atom -> channel_delta_cutoff; 
  double cutoff;


  for (n = 0; n<nlocal; n++){
    //    if (atype[n] != CHANNEL_ATOM_TYPE ) continue;
    if (atype[n] != ISOLATED_CHANNEL_ATOM_TYPE ) continue;
    channel_atom =n;
    rock_atom1 = atom->map(bond_rock_channel_atom[channel_atom][0]);
    rock_atom2 = atom->map(bond_rock_channel_atom[channel_atom][1]);
    
    delta_x = x[rock_atom1][0] - x[rock_atom2][0];
    delta_y = x[rock_atom1][1] - x[rock_atom2][1];
    delta_z = x[rock_atom1][2] - x[rock_atom2][2];
    dist = sqrt (delta_x *delta_x + delta_y *delta_y +delta_z *delta_z);

    max_bond_num = num_bond_channel_atom[channel_atom];    

    check = 0; // check if this channel is connected to wet channel atoms
    for (m = 0;m <max_bond_num;m++){
      temp_channel_atom = atom->map(bond_channel_atom[channel_atom][m]);
      check = (check || atype[temp_channel_atom] == CONNECTED_CHANNEL_ATOM_TYPE );
    }

    if (check == 1){
      // this dry channel atom becomes a wet channel atom now
      atype[channel_atom] = CONNECTED_CHANNEL_ATOM_TYPE;
      channel_width[channel_atom] = 0.0;
    } 

  }

}

/*==================================
(1) if spring extension is greater than cutoff, this spring breaks (SPRING_BOND_TYPE becomes DRY_CHANNEL_BOND_TYPE)
the corresponding channel atom becomes ISOLATED_CHANNEL_ATOM_TYPE
(2) in order to avoid a sudden channel walls movement as a dry crack becomes a wet one, we assume that the corresponding bond of a CONNECTED_CHANNEL_ATOM_TYPE bond is DRY_CHANNEL_BOND_TYPE by default. if this channel is filled with water (i.e., channel_width +1 >= dry crack width)
 this DRY_CHANNEL_BOND_TYPE becomes WET_CHANNEL_BOND_TYPE;
  ==================================*/

void FixLubricationUpdateBC::bond_break()
{
  int n, m, rock_atom1, rock_atom2;
  int nlocal = atom->nlocal;
  int *atype = atom -> type;
  int *num_bond = atom-> num_bond;
  int **bond_atom = atom-> bond_atom;
  double atom1x, atom1y, atom1z;
  double atom2x, atom2y, atom2z;
  double delta_x, delta_y, delta_z;
  double dist;
  double **x = atom->x;
  double **x0 = atom->x0;
  int **bond_type = atom -> bond_type;
  int channel_atom;
  int nstep = update->ntimestep;
  double *channel_width = atom->channel_width;
  int BC_cond; // unbreakable bonds outside this region
  double *channel_delta_cutoff = atom -> channel_delta_cutoff;
  double cutoff;

  for (n =0; n <nlocal; n++) {
    if (atype[n] != ROCK_ATOM_TYPE ) continue; // this is not a rock atom
    rock_atom1 = n;
    for (m = 0; m < num_bond[rock_atom1]; m++) {
      rock_atom2 = atom->map(bond_atom[rock_atom1][m]);
      if (atype[rock_atom2] != ROCK_ATOM_TYPE ) {fprintf(screen, "warning!! this should be a rock atom!!\n"); continue;}
      
          atom1x=x[rock_atom1][0];
	  atom1y=x[rock_atom1][1];
	  atom1z=x[rock_atom1][2];
	  
	  atom2x=x[rock_atom2][0];
	  atom2y=x[rock_atom2][1];
	  atom2z=x[rock_atom2][2];
	  
	  delta_x = atom1x - atom2x;
	  delta_y = atom1y - atom2y;
	  delta_z = atom1z - atom2z;
	  
	  dist = sqrt (delta_x*delta_x +delta_y*delta_y + delta_z*delta_z);
	  /*
	  if(dist > 1.1) {
	    fprintf(screen, "!!!!! rock takes apart, dist = %f too big !!!\n", dist );
	    error->one(FLERR,"rock takes apart!");
	    //	    error->done();	    //	    MPI_Finalize();	    //	    exit(1);
	  }
	  */

	  channel_atom = find_channel_atom(rock_atom1, rock_atom2);

	  cutoff = channel_delta_cutoff[channel_atom];

	  if (bond_type[rock_atom1][m]== SPRING_BOND_TYPE && dist >cutoff){

	    BC_cond = 1;
	    if ( (x0[channel_atom][0]- BC_xlo) * (x0[channel_atom][0]- BC_xhi) > 0 ) { BC_cond = 0;} // outside the region of fracturing system, BC_cond = 0 means unbreakable reagion.
	    if ( (x0[channel_atom][1]- BC_ylo) * (x0[channel_atom][1]- BC_yhi) > 0 ) { BC_cond = 0;} // outside the region of fracturing system, BC_cond = 0 means unbreakable reagion.
	    if ( (x0[channel_atom][2]- BC_zlo) * (x0[channel_atom][2]- BC_zhi) > 0 ) { BC_cond = 0;} // outside the region of fracturing system, BC_cond = 0 means unbreakable reagion.

	    if (BC_cond == 0 ) continue; // outside the main fracture system
	      
	    fprintf(screen, "ntimestep %d  bond break dist %f > cutoff %f at atom1xyz,atom2xyz (%.0f %.0f %.0f), (%.0f %.0f %.0f) \n", nstep, dist, cutoff, x0[rock_atom1][0], x0[rock_atom1][1],x0[rock_atom1][2], x0[rock_atom2][0], x0[rock_atom2][1],x0[rock_atom2][2]);

	    bond_type[rock_atom1][m] = DRY_CHANNEL_BOND_TYPE;
	    atype[channel_atom] = ISOLATED_CHANNEL_ATOM_TYPE;
	  }

	  // channel bond type upate 

	  if ( bond_type[rock_atom1][m] == DRY_CHANNEL_BOND_TYPE && atype[channel_atom] == CONNECTED_CHANNEL_ATOM_TYPE ){ //dry channel bond becomes wet channel bond if it is a wet channel atom
	    x[channel_atom][0] = 0.5*(x[rock_atom1][0] + x[rock_atom2][0] );
	    x[channel_atom][1] = 0.5*(x[rock_atom1][1] + x[rock_atom2][1] );
	    x[channel_atom][2] = 0.5*(x[rock_atom1][2] + x[rock_atom2][2] );
	    if (channel_width[channel_atom]+1 >dist) {
	      bond_type[rock_atom1][m] = WET_CHANNEL_BOND_TYPE;
	      fprintf(screen, "ntimestep %d  wet channel created, dist %f > cutoff %f at x y z %f %f %f \n", nstep, dist, cutoff, x0[channel_atom][0], x0[channel_atom][1], x0[channel_atom][2]);
	    }
	  }
  
    }
  
  }

}


/*------------------------------------------------------*/
int FixLubricationUpdateBC::pack_comm(int n, int *list, double *buf,
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

void FixLubricationUpdateBC::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){
    partner[i] = static_cast<int> (buf[m++]);
  }
}
/* ---------------------------------------------------------------------- */

int FixLubricationUpdateBC::pack_reverse_comm(int n, int first, double *buf)
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

void FixLubricationUpdateBC::unpack_reverse_comm(int n, int *list, double *buf)
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

double FixLubricationUpdateBC::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*sizeof(int);
  return bytes;
}
