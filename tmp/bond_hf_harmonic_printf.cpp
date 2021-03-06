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
#include "stdlib.h"
#include "bond_hf_harmonic.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "fix_lubrication_update_BC.h"
#include "types.h"
#include "region_block.h"
#include <iostream>
#include "update.h"


using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

BondHfHarmonic::BondHfHarmonic(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondHfHarmonic::~BondHfHarmonic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(k1);
    memory->destroy(beta);

  }
}

/* ---------------------------------------------------------------------- */

void BondHfHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **x0 = atom->x0;
  double **f = atom->f;
  double **f0 = atom->f0;
  double **v = atom->v;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double pot;

  int *tag = atom->tag;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;

  int atom1, atom2;
  double k_parallel, k_perp;
  double Delta[3],Delta_mag;
  double delij[3],vdiff[3],vmag;
  double dot_delij_Delta;
  double fijx,fijy,fijz;
  double fdamping_x, fdamping_y, fdamping_z;
  double fsum_x,fsum_y,fsum_z;
  int ntimestep = update->ntimestep;

  FILE *file = fopen ("test.txt", "w");


  //////////\\\\\\\\
    //  FOOX = 0.0;
    //FOOY = 0.0;
    //FOOZ = 0.0;
  ///////////\\\\\\\\

  for (n = 0; n<atom->nmax;n++){
    f0[n][0]= 0.0;
    f0[n][1]= 0.0;
    f0[n][2]= 0.0;
  }


  for (n = 0; n < nbondlist; n++) {
    atom1 = bondlist[n][0];
    atom2 = bondlist[n][1];
    type = bondlist[n][2];


    delx = x[atom1][0] - x[atom2][0]; // lammps' delta definitions
    dely = x[atom1][1] - x[atom2][1];
    delz = x[atom1][2] - x[atom2][2];
    
    k_parallel = k1[type];
    k_perp = k_parallel/3.0;
    
    
    Delta[0] = x0[atom2][0]-x0[atom1][0]; //The vector connecting neighbors originally
    Delta[1] = x0[atom2][1]-x0[atom1][1]; 
    Delta[2] = x0[atom2][2]-x0[atom1][2]; 
    Delta_mag = sqrt(Delta[0]*Delta[0]+Delta[1]*Delta[1]+Delta[2]*Delta[2]);
    
    /*-------  our delij definitions -----*/
    delij[0] = x[atom2][0]-x[atom1][0]-Delta[0];//Relative displacement of mass points
    delij[1] = x[atom2][1]-x[atom1][1]-Delta[1];
    delij[2] = x[atom2][2]-x[atom1][2]-Delta[2];
    
    rsq = (delij[0]+Delta[0])*(delij[0]+Delta[0])+(delij[1]+Delta[1])*(delij[1]+Delta[1])+(delij[2]+Delta[2])*(delij[2]+Delta[2]);
    r = sqrt(rsq);
    dr = r - r0[type];
    
    dot_delij_Delta=(delij[0]*Delta[0]+delij[1]*Delta[1]+delij[2]*Delta[2])/Delta_mag;
    
    
    
    // ------------- damping forces ---------- 
    // relative translational velocity
    //    vr1 = v[atom1][0] - v[atom2][0];
    //    vr2 = v[atom1][1] - v[atom2][1];
    //    vr3 = v[atom1][2] - v[atom2][2];
    // normal component
    //    rsqinv = 1.0/rsq;
    //    vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
    //    damp =  -1.0*beta[type]*vnnr;
    
    fijx = (k_parallel-k_perp)*dot_delij_Delta*Delta[0]/Delta_mag+k_perp*delij[0];//the contact force_x between i and j 
    fijy = (k_parallel-k_perp)*dot_delij_Delta*Delta[1]/Delta_mag+k_perp*delij[1];
    fijz =  (k_parallel-k_perp)*dot_delij_Delta*Delta[2]/Delta_mag+k_perp*delij[2];
    
    vdiff[0] = v[atom2][0]-v[atom1][0];
    vdiff[1] = v[atom2][1]-v[atom1][1];
    vdiff[2] = v[atom2][2]-v[atom1][2];
    vmag = (vdiff[0]*Delta[0]+vdiff[1]*Delta[1]+vdiff[2]*Delta[2])/Delta_mag;
    fdamping_x = beta[type]*vmag*Delta[0]/Delta_mag;
    fdamping_y = beta[type]*vmag*Delta[1]/Delta_mag;
    fdamping_z = beta[type]*vmag*Delta[2]/Delta_mag;
    // force & energy
    
    
    if (r <= 0.0){
      fsum_x = 0.0;
      fsum_y = 0.0;
      fsum_z = 0.0;
      continue;
    }
    
    if (type == SPRING_BOND_TYPE ){ // spring, restoring forces
      fsum_x = fijx + fdamping_x;
      fsum_y = fijy + fdamping_y;
      fsum_z = fijz + fdamping_z;
      /*      
      if ((update->ntimestep % 100 )==0 && (update->ntimestep > 10 )){
	//	fprintf(screen,"==fx,fy,fz %f %f %f fd_x,fd_y,fd_z %f %f %f, %d \n",1.e6*fijx,1.e6*fijy,1.e6*fijz,1.e6*fdamping_x,1.e6*fdamping_y,1.e6*fdamping_z);
	fprintf(screen,"==f_d/f[x],[y],[z] %f %f %f  \n",fdamping_x/(fijx+1e-6),fdamping_y/(fijy+1e-6),fdamping_z/(fijz+1e-6));
      }
      */
      //      fprintf(screen, "fdamping_x %10.2f beta %f vmag %f type %d\n",fdamping_x,beta[type],vmag,type);
    }
    else if (type == DRY_CHANNEL_BOND_TYPE ){ //|| type == WET_CHANNEL_BOND_TYPE ) {
      if (r > r0[type]){ // broken bond, zero forces under tensile.
	fsum_x = 0.0;
	fsum_y = 0.0;
	fsum_z = 0.0;
      }
      else{ //broke bond but contacted, so it behaves like granular particles under compression.
	fsum_x = fijx + fdamping_x;
	fsum_y = fijy + fdamping_y;
	fsum_z = fijz + fdamping_z;
      }
    }
    else if (type == WET_CHANNEL_BOND_TYPE ){
      fsum_x = 0.0;
      fsum_y = 0.0;
      fsum_z = 0.0;
    }
    else {fprintf(screen, "!!!!!!Warning!!!!!!\n");}

    /////////\
    //    std::cout << ntimestep << '\n';
    /////////\

    if (eflag) ebond =  k1[type]*dr*dr;
    
     if ( (x0[atom1][0]==global_xlo) || (x0[atom1][0]==global_xhi) || (x0[atom1][1]==global_ylo) || (x0[atom1][1]==global_yhi) || (x0[atom1][2]==global_zlo) || (x0[atom1][2]==global_zhi) ) {

      f[atom1][0] =0;
      f[atom1][1] =0;
      f[atom1][2] =0;
      
      f0[atom1][0] =0;
      f0[atom1][1] =0;
      f0[atom1][2] =0;
      ///////////\\\\\\\\\\\\\\\
      if ( (x0[atom1][0]==9) && (x0[atom1][1]==102) && (x0[atom1][2]==4) && (x0[atom2][0]==9) && (x0[atom2][1]==102) && (x0[atom2][2]==3) ){
	//   	      FOOX = fsum_x;
	//	      FOOY = fsum_y;
	//	      FOOZ = fsum_z;
	/*	FOOX = x0[atom1][0];
	FOOY = x0[atom1][1];
	FOOZ = x0[atom1][2];
	*/
	if (ntimestep % 500 == 0) {
	  //	  std::cout << x0[atom1][0] <<" "<< x0[atom1][1] << " " << x0[atom1][2] << std::endl;
	  fprintf(file, "%d %g %g %g \n",ntimestep, fsum_x, fsum_y, fsum_z);
	} 

	////////
	
	
      } // end of if of (x0[atom1][0]==9)
      //////////\\\\\\\\\\\\\\\\
    }
    
    else if (newton_bond || atom1 < nlocal) {
      
      f[atom1][0] += fsum_x;
      f[atom1][1] += fsum_y;
      f[atom1][2] += fsum_z;
      
      f0[atom1][0] += fsum_x;
      f0[atom1][1] += fsum_y;
      f0[atom1][2] += fsum_z;

    }
    
     if ( (x0[atom2][0]==global_xlo) || (x0[atom2][0]==global_xhi) || (x0[atom2][1]==global_ylo) || (x0[atom2][1]==global_yhi) || (x0[atom2][2]==global_zlo) || (x0[atom2][2]==global_zhi) ) {
    
      f[atom2][0] =0;
      f[atom2][1] =0;
      f[atom2][2] =0;
      
      f0[atom2][0] =0;
      f0[atom2][1] =0;
      f0[atom2][2] =0;
      // ///////\\\\\\\\\\\\\\\\      
      if ( (x0[atom1][0]==9) && (x0[atom1][1]==102) && (x0[atom1][2]==3) && (x0[atom2][0]==9) && (x0[atom2][1]==102) && (x0[atom2][2]==4) ){
      // 	//   	      FOOX = -fsum_x;
      // 	//	      FOOY = -fsum_y;
      // 	//	      FOOZ = -fsum_z; 
	//	FOOX = x0[atom2][0];
	//FOOY = x0[atom2][1];
	//FOOZ = x0[atom2][2];
	if (ntimestep % 500 == 0) {
	  //	  std::cout << x0[atom2][0] <<" "<< x0[atom2][1] << " " << x0[atom2][2] << std::endl;
	  fprintf(file, "%d %g %g %g \n",ntimestep, fsum_x, fsum_y, fsum_z);
	}
	
      }
      //////\\\\\\\\\\\\\\\\\\\
    }

    else if (newton_bond || atom2 < nlocal) {
      
      f[atom2][0] -= fsum_x;
      f[atom2][1] -= fsum_y;
      f[atom2][2] -= fsum_z;
      
      f0[atom2][0] -= fsum_x;
      f0[atom2][1] -= fsum_y;
      f0[atom2][2] -= fsum_z;

    }
     ///    
     fclose(file);
    
    if (evflag) ev_tally(atom1,atom2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
  
  rock_force_update();
  comm->forward_comm();
  comm->forward_comm();
}

/* ---------------------------------------------------------------------- */

void BondHfHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;
  memory->create(r0,n+1,"bond:r0");
  memory->create(k1,n+1,"bond:k1");
  memory->create(beta,n+1,"bond:beta");


  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHfHarmonic::coeff(int narg, char **arg)
{
  if (narg !=4) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double r0_one = force->numeric(FLERR,arg[1]);
  double k1_one = force->numeric(FLERR,arg[2]);
  double beta_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    r0[i] = r0_one;
    k1[i] = k1_one;
    beta[i] = beta_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondHfHarmonic::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHfHarmonic::write_restart(FILE *fp)
{
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k1[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&beta[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHfHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&k1[1],sizeof(double),atom->nbondtypes,fp);
    fread(&beta[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k1[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&beta[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondHfHarmonic::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k1[type] * dr;
  double pot;
  fforce = 0;
  if (r > 0.0) fforce = - k1[type]* (r-r0[type]);
  pot = 0.5* k1[type]*dr*dr;
  return pot;
}


/*------------------------------------------------------
  update channel pressure info
------------------------------------------------------*/

void BondHfHarmonic::rock_force_update(){
  int nlocal = atom->nlocal;
  double **x0 = atom->x0;
  double **v = atom->v;
  int *tag = atom->tag;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int newton_bond = force->newton_bond;
  int **bond_type = atom->bond_type;
  double *channel_pressure = atom->channel_pressure;

  int i,atom1, atom2;
  double atom1x, atom1y,atom1z, atom2x, atom2y,atom2z;
  double channelx, channely,channelz;
  double delta1x,delta1y,delta1z, deltamag;
  double delta2x,delta2y,delta2z;
  double f_atom1_x, f_atom1_y, f_atom1_z;
  double f_atom2_x, f_atom2_y, f_atom2_z;
  double **f = atom->f;
  int type;

  for (atom1 = 0; atom1 < nlocal; atom1++) {
    for (i = 0; i < num_bond[atom1]; i++) {
      atom2 = atom->map(bond_atom[atom1][i]);
      type = bond_type[atom1][i];
      if (type != WET_CHANNEL_BOND_TYPE ) continue;// this is not a wet channel
      if (newton_bond || atom1< atom2){
	
	atom1x=x0[atom1][0];
	atom1y=x0[atom1][1];
	atom1z=x0[atom1][2];

	atom2x=x0[atom2][0];
	atom2y=x0[atom2][1];
	atom2z=x0[atom2][2];
	
	channelx = 0.5* (atom1x + atom2x);
	channely = 0.5* (atom1y + atom2y);
	channelz = 0.5* (atom1z + atom2z);
	
	
	delta1x = atom1x - channelx;
	delta1y = atom1y - channely;
	delta1z = atom1z - channelz;
	
	delta2x = -delta1x;
	delta2y = -delta1y;
	delta2z = -delta1z;
	deltamag= sqrt( delta1x*delta1x + delta1y*delta1y + delta1z*delta1z );
	
	f_atom1_x = f[atom1][0];
	f_atom1_y = f[atom1][1];
	f_atom1_z = f[atom1][2];
	
	
	f_atom2_x = f[atom2][0];
	f_atom2_y = f[atom2][1];
	f_atom2_z = f[atom2][2];

	
	if (fabs(delta1x) > 1.e-6) {
	  if (newton_bond || atom1 < nlocal) {
	    f[atom1][0] = 0.0;
	    v[atom1][0] = 0.0;
	  }
	  if (newton_bond || atom2 < nlocal) {
	    f[atom2][0] = 0.0;
	    v[atom2][0] = 0.0;
	  }
	  
	}
	else if (fabs(delta1y) > 1.0e-6){
	  if (newton_bond || atom1 < nlocal) {
	    f[atom1][1] = 0.0;
	    v[atom1][1] = 0.0;
	  }
	  if (newton_bond || atom2 < nlocal) {
	    f[atom2][1] = 0.0;
	    v[atom2][1] = 0.0;
	  }
	  //	fprintf(screen,"set force zero \n");
	}
	else if (fabs(delta1z) > 1.0e-6){
	  if (newton_bond || atom1 < nlocal) {
	    f[atom1][2] = 0.0;
	    v[atom1][2] = 0.0;
	  }
	  if (newton_bond || atom2 < nlocal) {
	    f[atom2][2] = 0.0;
	    v[atom2][2] = 0.0;
	  }
	  
	}
      }
    }
  }
}

