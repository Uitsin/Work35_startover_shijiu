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

#include "string.h"
#include "math.h"
#include "fix_lubrication_update_PKN.h"
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

FixLubricationUpdatePKN::FixLubricationUpdatePKN(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 12) error->all(FLERR,"Illegal fix lubrication/update command"); // ID group-ID fixchannelini
  //  cutoff = atof(arg[3]);
  injection_x = atof(arg[3]);
  injection_y = atof(arg[4]);
  injection_z = atof(arg[5]);
  PKN_X =  atof(arg[7]);
  PKN_Y =  atof(arg[9]);
  PKN_Z =  atof(arg[11]);

  partner = NULL;
  nmax =0;
}



FixLubricationUpdatePKN::~FixLubricationUpdatePKN()
{
  memory->destroy(partner);
}


/*============
  setmask
  ============*/
int FixLubricationUpdatePKN::setmask()
{
  int mask = 0;
   mask |= PRE_FORCE;
   mask |= FINAL_INTEGRATE;
  return mask;
}


/*==================================
  called before force routine
  ==================================*/

//void FixLubricationUpdatePKN::pre_force(int vflag)
void FixLubricationUpdatePKN::final_integrate()

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

  lubrication();
  comm->forward_comm();


  neighbor->build_topology();
  channel_update();
  comm->forward_comm();



}


/*==================================
  calculate channel pressure
  ==================================*/
void FixLubricationUpdatePKN::cal_channel_pressure()
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

    if (pressure < 0) pressure = 0; // this channel is partially filled with water

    channel_pressure[channel_atom] = pressure;

    if ((update->ntimestep %50000)==0) fprintf(screen, "channel pressure is %.2f MPa at x y z %f %f %f  \n",channel_pressure[channel_atom]*1.e-6,x0[channel_atom][0],x0[channel_atom][1],x0[channel_atom][2]);
    
  }

}


/*======================================
  routine used to check channel pressure only 
  ======================================*/
/*
void FixLubricationUpdatePKN:: check_channel_pressure()
{
  int n;
  int nbondlist = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;
  int rock_atom1,rock_atom2,btype;
  int channel_atom;
  int *atype = atom->type;
  double *channel_pressure = atom->channel_pressure;
  double **x0 = atom->x0;
  for (n = 0; n < nbondlist; n++) {
    rock_atom1 = bondlist[n][0];
    rock_atom2 = bondlist[n][1];
    btype = bondlist[n][2];
    if (btype != WET_CHANNEL_BOND_TYPE ) continue;
    channel_atom = find_channel_atom(rock_atom1, rock_atom2);
    //    fprintf(screen, "===== check atype of this channel is %d (3?), channel_atom= %d at x y %f %f,  pressure %f \n",atype[channel_atom],channel_atom, x0[channel_atom][0],x0[channel_atom][1],channel_pressure[channel_atom]);
    
  }
}
*/



/*========================================
  channel flow based on lubrication approximation
  ========================================*/

void FixLubricationUpdatePKN::lubrication(){

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
  int rock_atom1, rock_atom2;
  int nstep = update->ntimestep;
  double ntime = update->dt;
  double dw;

  on_twelve_mu_L =1.0/(12.0*mu*L);// need to be rescaled if lattice spacing is not 1 meter.
  for (n=0; n<nlocal;n++){
    if (atype[n] != CONNECTED_CHANNEL_ATOM_TYPE) continue;
    channel_width_old[n]= channel_width[n];
  }

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

      if (channel_pi > channel_pj){local_w = channel_wi;}
      else {local_w = channel_wj;}
      vij = (channel_pj-channel_pi)*local_w*local_w*local_w*on_twelve_mu_L;
      //      channel_width[channel_atomi]+=vij *update->dt;
      dw = vij * update->dt;
      channel_width[channel_atomi]+=0.5*dw;
      channel_width[channel_atomj]-=0.5*dw;
      //channel_width[channel_atomi]+=dw;
    }
  }
}


/*========================================
  channel atom update based on new channel width
  ========================================*/
void FixLubricationUpdatePKN::channel_update(){

  int n;
  int nbondlist = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;
  int rock_atom1, rock_atom2,channel_atom;
  int btype;
  int *atype = atom->type;
  double  *channel_width = atom->channel_width;
  double channel_w;
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

int FixLubricationUpdatePKN::find_channel_atom(int rock_atom1, int rock_atom2){
  
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






void FixLubricationUpdatePKN::new_channel()
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
  int connection_check[6];//={0,0,0,0,0,0};
  
  int check;
  int m;
  int PKN_cond;
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


    /*
    cutoff = channel_delta_cutoff[channel_atom];

    //    if (dist >cutoff && check ==1 ){
    if (dist >cutoff){
      
      PKN_cond = 1 ;
           
      if ( fabs(x0[channel_atom][0]- injection_x)> 0.5 *PKN_X )  {
	PKN_cond = 0; 
      }
      if ( fabs(x0[channel_atom][1]- injection_y)> 0.5 *PKN_Y )  {
	PKN_cond = 0; 
      }
      
      if ( fabs(x0[channel_atom][2]- injection_z)> 0.5 *PKN_Z )  {
	PKN_cond = 0; 
      }
      
      if (PKN_cond == 0 ) continue; // outside the main fracture system

      if (check == 1) { // this atom becomes a wet channel atom now
	atype[channel_atom] = CONNECTED_CHANNEL_ATOM_TYPE;
	channel_width[channel_atom] = 0.0;
      } 
      //      else if (check == 0) {atype[channel_atom] = ISOLATED_CHANNEL_ATOM_TYPE;} // this atom becomes a dry channel atom now

    }
    */

  }

}




void FixLubricationUpdatePKN::bond_break()
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
  int PKN_cond; // add restrictions to fix fracture height   
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
	  if(dist > 1.1) {
	    fprintf(screen, "!!!!! rock takes apart, dist = %f too big !!!\n", dist );
	    error->one(FLERR,"rock takes apart!");
	    //	    error->done();	    //	    MPI_Finalize();	    //	    exit(1);
	  }


	  channel_atom = find_channel_atom(rock_atom1, rock_atom2);

	  //	  fprintf(screen, "injecting at x y z %.1f %.1f %.1f PKN X Y Z  %.1f %.1f %.1f\n",injection_x, injection_y, injection_z, PKN_X, PKN_Y, PKN_Z);
	  //	  fprintf(screen, "fixed pkn_y  is %f\n ", PKN_Y);
	  cutoff = channel_delta_cutoff[channel_atom];

	  if (bond_type[rock_atom1][m]== SPRING_BOND_TYPE && dist >cutoff){
	    PKN_cond = 1 ;

	    
	    if ( fabs(x0[channel_atom][0]- injection_x)> 0.5 *PKN_X )  {
	      PKN_cond = 0; 
	       if ((update->ntimestep %50000)==0) fprintf(screen, "PKN x-constraint is on at x y z %.0f %.0f %.0f \n", x0[channel_atom][0],x0[channel_atom][1], x0[channel_atom][2]);
	    }
	    if ( fabs(x0[channel_atom][1]- injection_y)> 0.5 *PKN_Y )  {
	      PKN_cond = 0; 
	       if ((update->ntimestep %50000)==0) fprintf(screen, "PKN y-constraint is on at x y z %.0f %.0f %.0f \n", x0[channel_atom][0],x0[channel_atom][1], x0[channel_atom][2]);
	    }

	    if ( fabs(x0[channel_atom][2]- injection_z)> 0.5 *PKN_Z )  {
	      PKN_cond = 0; 
	       if ((update->ntimestep %50000)==0) fprintf(screen, "PKN z-constraint is on at x y z %.0f %.0f %.0f \n", x0[channel_atom][0],x0[channel_atom][1], x0[channel_atom][2]);
	    }


	    if (PKN_cond == 0 ) continue; // outside the main fracture system
	      
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
int FixLubricationUpdatePKN::pack_comm(int n, int *list, double *buf,
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

void FixLubricationUpdatePKN::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++){
    partner[i] = static_cast<int> (buf[m++]);
  }
}
/* ---------------------------------------------------------------------- */

int FixLubricationUpdatePKN::pack_reverse_comm(int n, int first, double *buf)
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

void FixLubricationUpdatePKN::unpack_reverse_comm(int n, int *list, double *buf)
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

double FixLubricationUpdatePKN::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*sizeof(int);
  return bytes;
}
