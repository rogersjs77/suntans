/*
 * File: sources.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Right-hand sides for momentum and scalars, in the form rhs = dt*SOURCE at time step n.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "phys.h"
#include "sources.h"
#include "memory.h"

void MomentumSource(REAL **usource, gridT *grid, physT *phys, propT *prop, averageT *average) {

  int j, jptr, nc1, nc2, k;
  REAL Coriolis_f, ubar, depth_face;

/* SCS realistic example flows
  REAL tauL = 14*12.42*3600, ULm = 0.10, ULz=0; 
  REAL TM2 = 12.42*3600, UHtide = 0.35, Phitide = 1, UHiw=0.50, Phiiw = 0, tauW=0.25;
  REAL alpha1 = -0.7198, alpha2=-1.6810, delta=300, drho=5.5, h2=2500;
  REAL F_hatx=1, F_haty=0; 
*/
  REAL ULm = prop->ULm, ULz=prop->ULz, UHtide = prop->UHtide, UHiw = prop->Uiw;
  REAL VLm = prop->VLm, VLz=prop->VLz;
  REAL TauL = prop->TauL*3600, TM2 = prop->TM2*3600; 
  REAL Phitide = prop->Phitide, Phiiw = prop->Phiiw; 
  REAL alpha1 = prop->alpha1, alpha2=prop->alpha2, delta=prop->delta, drho=prop->drho, h2=prop->h2;
  REAL D=prop->D;
  REAL xe, ye, z, fz, fzb;
  REAL D_hat, F_hat, k_hat;
  REAL u_bar, u_tilde, U_bar, U_tilde, u_bar_prime, u_tilde_prime;
  REAL u_LS_bar, U_LS_tilde, u_LS_tilde_prime, v_LS_bar, V_LS_tilde, v_LS_tilde_prime;
  REAL tr_hat = 1-exp(-prop->rtime/prop->thetaramptime);// rampup factor
  REAL deltaD = prop->sponge_distance, tauD = prop->sponge_decay;
  
  // compute iw wavenumbers
  REAL kiw_x = 2*PI/prop->lambda*prop->Fhat_x;
  REAL kiw_y = 2*PI/prop->lambda*prop->Fhat_y;
  REAL kiw = sqrt(pow(kiw_x,2)+pow(kiw_y,2));

  // compute tide wavenumbers
  REAL Ctide = sqrt(D*9.81);
  REAL ktide = 2*PI/(TM2*Ctide);
  REAL ktide_x = ktide*prop->Fhat_x;
  REAL ktide_y = ktide*prop->Fhat_y;

  //printf("sources xmin=%f, xmax=%f\n", xmin,xmax);
  //ComputeAlphaw(grid, prop, average);
  // go through each edge 
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {

	  // get initial parameters
    j = grid->edgep[jptr];     
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    xe = 0.5*(grid->xv[nc1]+grid->xv[nc2]); //edge distance
    ye = 0.5*(grid->yv[nc1]+grid->yv[nc2]); //edge distance
   
    // wave flux direction
    if( (prop->Fhat_x*rn1[j]+prop->Fhat_y*rn2[j]) < 0){   // Flux dot rnormal is negative = incoming wave
      F_hat=1;
    }else{
      F_hat=0;
    }

    // wave velocity direction
    // REAL k_hatx=1, k_haty=0;
    // if( (k_hatx*rn1[j]+k_haty*rn2[j]) < 0){   // Flux dot rnormal is negative = incoming wave
    //   k_hat=1;
    // }else{
    //   k_hat=0;
    // }

    // wavemaker location, radial from edge of sponge layer plus dW width
    // tanh function wavemaker with width dW, side rolloff=dWside
    // REAL dWtop=dW-2*dWside;
    // REAL rWmid=prop->sponge_distance+dW/2;
    // REAL rW1=rWmid-dWtop/2-dWside/2;
    // REAL rW2=rWmid+dWtop/2+dWside/2;
    // W_hat = 0.5*tanh(2*(rSponge[j]-rW1)/dWside)+0.5*tanh(-2*(rSponge[j]-rW2)/dWside);
    // // exponental wavemaker
    // REAL a1=0.5, lambda=160e3;
    // W_hat = exp(-8*pow(rSponge[j]-rWmid,2)/pow(a1*lambda,2));
    // rectangular wavemaker
    // if(prop->sponge_distance<rSponge[j] && (prop->sponge_distance+dW)>rSponge[j]){
    //   W_hat=1;
    // }else{
    //   W_hat=0;
    // }

    // sponge layer damping coeff, only from boundary to deltaD
    D_hat=0;
    if(rSponge[j]<deltaD){
      D_hat = exp(-4.0*rSponge[j]/deltaD);
    }    

    /* velocity decomposition */
    //u = u_bar + u_tilde (separate in time, low pass, high pass)
    //u = U + u_prime (separate in space, baroclinic, barotropic)
    //u = (U_bar + u_bar_prime) + (U_tilde + u_tilde_prime); (bt mean, bc mean, bt tide, internal waves)
    //u_bar = average->u_avg;
    //u_tilde = phys->u - average->u_avg;
    //U_bar = depth_avg(u_bar);
    //U_tilde = depth_avg(u_tilde;
    //u_bar_prime = u_bar - U_bar;
    //u_tilde_prime = u_tilde - U_tilde;

    /* compute depth avg velocity */ 
    U_bar = U_tilde = 0;
    depth_face = 0;
    fzb=0;
    z=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      z-=grid->dz[k]/2; // z on faces is at midpoint
      u_bar = average->u_avg[j][k];
      u_tilde = phys->u[j][k] - u_bar;
      U_bar += grid->dzf[j][k]*u_bar;
      U_tilde += grid->dzf[j][k]*u_tilde;
      depth_face += grid->dzf[j][k];
      z-=grid->dz[k]/2;
    }
    U_bar/=depth_face; // mean barotropic
    U_tilde/=depth_face; // mean baroclinic
    
    /* Apply forcing, cycle through each z layer*/
    z=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      z-=grid->dz[k]/2; // z on faces is at midpoint

      // compute velocity decomposition (edge normal)
      // u = (U_bar + u_tilde_bar) + (U_prime + u_tilde_prime);
      u_bar = average->u_avg[j][k];
      u_tilde = phys->u[j][k] - u_bar;
      u_bar_prime = u_bar - U_bar;
      u_tilde_prime = u_tilde - U_tilde;

      // compute forcing velocities (large scale)(u_LS,v_LS)
      // u_LS = u_LS_bar + u_LS_tilde (time decomposition)
      // fz = alpha1 * tanh((z+D-h2)/delta) + alpha2;
      fz = alpha1 * exp(-(-z+h2)/delta) + alpha2;
      u_LS_bar = ULm + ULz*fz;
      U_LS_tilde =  UHtide * (ktide_x/ktide) * cos(ktide_x*xe + ktide_y*ye - 2*PI*prop->rtime/TM2 + Phitide);
      u_LS_tilde_prime =  UHiw * fz * (kiw_x/kiw) * cos(kiw_x*xe + kiw_y*ye - 2*PI*prop->rtime/TM2 + Phiiw);
      
      v_LS_bar = VLm + VLz*fz;
      V_LS_tilde =  UHtide * (ktide_y/ktide) * cos(ktide_x*xe + ktide_y*ye - 2*PI*prop->rtime/TM2 + Phitide);
      v_LS_tilde_prime =  UHiw * fz * (kiw_y/kiw) * cos(kiw_x*xe + kiw_y*ye - 2*PI*prop->rtime/TM2 + Phiiw);
      
      /* This is the sponge layer */
      if(prop->sponge_distance) { 
        usource[j][k]-= prop->dt * (u_tilde_prime - 
          F_hat*tr_hat*u_LS_tilde_prime*grid->n1[j] - 
          F_hat*tr_hat*v_LS_tilde_prime*grid->n2[j])*D_hat/tauD;
      }
    
      /* this is the wave forcing term*/
      //usource[j][k]-= prop->dt*average->alphaw*du_wave_dt1*tr_hat*W_hat*F_hat;
        
      // high frequency forcing in middle for testing
      // if(xe>=Lmid1 && xe<=Lmid2) {
      //   usource[j][k]-= prop->dt*du_wave_dt1*tr_hat;
      // }

      /*Low frequency forcing over entire domain*/
      usource[j][k]-= prop->dt*(u_bar - tr_hat*u_LS_bar*grid->n1[j] - tr_hat*v_LS_bar*grid->n2[j])
              /TauL;

      /* extra forcing for periodic in y, geostrophic balance*/
      usource[j][k]-= -prop->dt*(u_LS_bar*grid->n2[j]*prop->Coriolis_f)*
               tr_hat;

      /* extra forcing for periodic in x, geostrophic balance*/
      // usource[j][k]-= prop->dt*(v_LS_bar*grid->n1[j]*prop->Coriolis_f)*
      //         tr_hat;          

      z-=grid->dz[k]/2;
    }
    // if(fabs(prop->rtime)>100000)
    // printf("j=%d, u_bar=%f, uLS*n1[j]=%f, vLS*n2[j]=%f, phys->u[j][0]=%f, usource[j][0]=%f \n", j, u_bar,u_LS_bar*grid->n1[j], v_LS_bar*grid->n2[j],phys->u[j][0],usource[j][0] );
  }

  /* Coriolis for a 2d problem */    

  if(prop->n==prop->nstart+1) {
    //printf("Initializing v Coriolis for a 2.5D simulation\n");
    v_coriolis = (REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"MomentumSource");
    for(j=0;j<grid->Ne;j++) {
      v_coriolis[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"MomentumSource");

      for(k=0;k<grid->Nke[j];k++)
        v_coriolis[j][k] = 0.0;
    }
  }

  // Hard-code coriolis here so that it can be zero in the main code
  Coriolis_f=0;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 
      
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
      
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      usource[j][k]+=prop->dt*Coriolis_f*(v_coriolis[j][k]*grid->n1[j]-
              InterpToFace(j,k,phys->uc,phys->u,grid)*grid->n2[j]);
  }

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 
      
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
       v_coriolis[j][k]-=prop->dt*Coriolis_f*InterpToFace(j,k,phys->uc,phys->u,grid);
  }

}

void KurtSource(gridT *grid, physT *phys, propT *prop, MPI_Comm comm) {
/* Kurt style forcing to maintain constant volume averaged velocity (and reach steady state fast)
Steps:
1) compute volume averaged velocity from u_hat, the velocity currently stored in phys->u[j][k]
2) compute forcing necessary to make the volume average of u_n+1 = u0
3) add forcing times dt to usource[j][k]
*/
int i, iptr, j, jptr, nc1, nc2, k; 
REAL uvol = 0, cellvol = 0, cumvol = 0, myuvol=0, mycumvol=0, u0=0, S;
// for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
//   j = grid->edgep[jptr];
//   nc1 = grid->grad[2*j];
//   nc1 = grid->grad[2*j+1];
//   cellarea = 0.5*(grid->Ac[nc1]+grid->Ac[nc2]);
//   depth_face = 0;
//   for(k=grid->etop[j];k<grid->Nke[j];k++) {
//     uvol+= phys->u[j][k]*grid->dzf[j][k];
//     depth_face += grid->dzf[j][k];
//   }
//   cumvol+=depth_face;
// }

  if(fabs(u0)>0){ // only run routine if specified velocity
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      for(k=grid->ctop[i];k<grid->Nk[i];k++){
          cellvol = grid->Ac[i]*grid->dzz[i][k];
          myuvol+=phys->uc[i][k]*cellvol;
          mycumvol+=cellvol;
          // cellvol=0; //not necessary
      }
    }

     MPI_Reduce(&mycumvol,&(cumvol),1,MPI_DOUBLE,MPI_SUM,0,comm);
     MPI_Bcast(&cumvol,1,MPI_DOUBLE,0,comm);
     MPI_Reduce(&myuvol,&(uvol),1,MPI_DOUBLE,MPI_SUM,0,comm);
     MPI_Bcast(&uvol,1,MPI_DOUBLE,0,comm);
     
    uvol/=cumvol;
    // rather than assigning each cell to receive uniform forcing (u0-uvol), we could account for hill's presence by
    // assigning each vertical column of cells to receive unifom depth integrated forcing ((u0-uvol)*D0/D(x))... but not doing that above

    S = (u0-uvol)/prop->dt;
    //note that this does not account for variable depth... could instead calculate column average velocity and do a S on that

    //write S to a new file if this is the first timestep
    if(prop->n==1+prop->nstart) {
      FILE *Sfp = fopen( "KurtS.txt" , "w" );
      fprintf (Sfp,"%e \n", S);
      fclose(Sfp);
    }
    //append S to file if n/ntout is an integer
    if(!(prop->n%prop->ntout)) {
      FILE *Sfp = fopen( "KurtS.txt" , "a" );
      fprintf (Sfp,"%e \n", S);
      fclose(Sfp);
    }

    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++){
      j = grid->edgep[jptr];

      // nc1 = grid->grad[2*j];
      // nc2 = grid->grad[2*j+1];
      // xe = 0.5*(grid->xv[nc1]+grid->xv[nc2]);
      // ye = grid->ye[j];
      // S = (u0-uvol)*D/ReturnDepth(xe,ye)/prop->dt;

      // update u with S/dt
      for(k=grid->etop[j];k<grid->Nke[j];k++){
        phys->u[j][k]+=(prop->dt*S)*grid->n1[j];
      }

    }
  }

}


/*
 * Function: HeatSource
 * Usage: HeatSource(grid,phys,prop,A,B);
 * --------------------------------------
 * Source terms for heat equation of the form
 *
 * dT/dt + u dot grad T = d/dz ( kappaT dT/dz) + A + B*T
 *
 * Assuming adiabatic top and bottom.  Horizontal advection is
 * explicit while all other terms use the theta method.
 *
 * Note that they must be set to zero if there is no heat
 * flux since these are temporary variables with values that may
 * be assigned elsewhere.
 *
 */
void HeatSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, metT *met, int myproc, MPI_Comm comm) {
  int i, k;
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      A[i][k]=B[i][k]=0;
}

/*
* Funtion: SaltSource
*/

void SaltSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, metT *met) {
  int i, k;
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      A[i][k]=B[i][k]=0;
}


/* ** not used anymore **
 * Function: ComputeAlphaw(grid,average,prop);
 * -------------------------------------
 * Computes the term alphaw which is the transfer
 * function coefficient for wave forcing
 *
 */
/*void ComputeAlphaw(gridT *grid, propT *prop, averageT *average, int myproc, int numprocs, MPI_Comm comm) {

  int j, jptr, nc1, nc2, k, proc;
  REAL sum_a, mysum_a, sum_kk, mysum_kk;
  MPI_Status status;
  REAL TauL = prop->TauL*3600, ULm = prop->ULm, ULz=prop->ULz; 
  REAL TM2 = prop->TM2*3600, UHtide = prop->UHtide, Phitide = prop->Phitide; 
  REAL UHiw = prop->Uiw, Phiiw = prop->Phiiw; 
  REAL alpha1 = prop->alpha1, alpha2=prop->alpha2, delta=prop->delta, drho=prop->drho, h2=prop->h2;
  REAL D=prop->D, F_hatx=prop->Fhat_x, F_haty=prop->Fhat_y;
  REAL dW=prop->dW;
  REAL xe, z, fz, fzb;
  REAL W_hat, D_hat, F_hat;
  REAL tr_hat = 1-exp(-prop->rtime/prop->thetaramptime);// rampup factor
  REAL deltaD = prop->sponge_distance, tauD = prop->sponge_decay;
  REAL u_rms_LS, u_rms, alphaw, kk, z_alphaw=prop->z_alphaw;

  //printf("start average->alphaw=%f\n", average->alphaw);

  // for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++){
  //   sum+=1;
  //   }
  //   printf("For loop, myproc=%d, %f\n",myproc,sum);


  // compute wave transfer coefficient alphaw
  // alphaw=0;
  // kk=0;

  // get z first (same for all cells)
  jptr=grid->edgedist[0];
  // get initial parameters
  j = grid->edgep[jptr];     
  // go to just below depth z_alphaw
  z=0;
  for(k=grid->etop[j];k<grid->Nke[j];k++) {
    z-=grid->dz[k]/2; // z on faces is at midpoint
      if( z < z_alphaw){
        break; // go to next face
      }
      z-=grid->dz[k]/2; 
    }  
    // printf("z=%f\n",z );


  // find relevant cells in each processor and add up sum_a and sum_kk
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {

    // get initial parameters
    j = grid->edgep[jptr];     
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    xe = 0.5*(grid->xv[nc1]+grid->xv[nc2]); //edge distance

    fz = alpha1 * tanh((z+D-h2)/delta) + alpha2;
    //printf("fz=%f, alpha1=%f, z=%f, D=%f, h2=%f, delta=%f, alpha2%f\n",fz,alpha1,z,D,h2,delta,alpha2);
    // theoretical rms of wave velocity, always positive
    u_rms_LS = fabs(UHiw*fz*grid->n1[j])*tr_hat;
    //u_rms_LS = fabs(UHiw*fz*grid->n1[j]) + fabs(0*fz*grid->n2[j]);
    // computed rms of wave velocity from variance
    u_rms = sqrt(2*average->uw_var[j][k]);

    // set constraints and compute alphaw only if cond met:
    // 1. location is inner 1/4 of wavemaker
    // 2. wave direction is not perpendicular to face
    // 3. wave flux is into domain
    if (rSponge[j]<(deltaD+1.5*dW) &&
      rSponge[j]>(deltaD+dW) &&
      u_rms>0.01 &&
      (F_hatx*rn1[j]+F_haty*rn2[j]) < 0){
      sum_a+=u_rms/u_rms_LS;
      //sum_a+=tanh(-u_rms/u_rms_LS+1)+1;
      sum_kk+=1;
      //printf("myproc=%d, xe=%f, z=%f, urms=%f, urms_LS=%f, fz=%f, sum_a=%f, sum_kk=%f\n",myproc,xe,z,u_rms,u_rms_LS,fz,sum_a,sum_kk);
    }
  }


  //write alphaw to a new file if this is the first timestep
  // if(prop->n==1+prop->nstart) {
  //   FILE *Sfp = fopen( "./data/alphaw.dat" , "w" );
  //   fprintf (Sfp,"%e %e \n", average->alphaw, prop->rtime);
  //   fclose(Sfp);
  // }else{
  // //append alphaw to file if n/ntout is an integer
  // // if(!(prop->n%prop->ntout)) {
  //   FILE *Sfp = fopen( "./data/alphaw.dat" , "a" );
  //   fprintf (Sfp,"%e %e \n", average->alphaw, prop->rtime);
  //   fclose(Sfp);
  // }


// send and receive sum_a
  if(myproc==0){
    for(proc=1;proc<numprocs;proc++){
      MPI_Recv(&mysum_a,1,MPI_DOUBLE,MPI_ANY_SOURCE,1,comm,&status);
      sum_a+=mysum_a;
      // printf("MPI_Recv a myproc=%d, proc=%d, mysum_a=%f, sum_a=%f, mysum_kk=%f, sum_kk=%f\n",myproc,proc,mysum_a,sum_a,mysum_kk,sum_kk);
    }
  }else{
    MPI_Send(&sum_a,1,MPI_DOUBLE,0,1,comm);
    // printf("MPI_Send a myproc=%d, mysum_a=%f, sum_a=%f, mysum_kk=%f, sum_kk=%f\n",myproc,mysum_a,sum_a,mysum_kk,sum_kk);
  }
  MPI_Bcast(&sum_a,1,MPI_DOUBLE,0,comm);
  // printf("MPI_Bcast a myproc=%d, sum_a=%f, sum_kk=%f\n",myproc,sum_a,sum_kk);

  // send receive sum_kk
  if(myproc==0){
    for(proc=1;proc<numprocs;proc++){
      MPI_Recv(&mysum_kk,1,MPI_DOUBLE,MPI_ANY_SOURCE,1,comm,&status);
      sum_kk+=mysum_kk;
      // printf("MPI_Recv kk myproc=%d, proc=%d, mysum_a=%f, sum_a=%f, mysum_kk=%f, sum_kk=%f\n",myproc,proc,mysum_a,sum_a,mysum_kk,sum_kk);
    }
  }else{
    MPI_Send(&sum_kk,1,MPI_DOUBLE,0,1,comm);
    // printf("MPI_Send kk myproc=%d, mysum_a=%f, sum_a=%f, mysum_kk=%f, sum_kk=%f\n",myproc,mysum_a,sum_a,mysum_kk,sum_kk);
  }
  MPI_Bcast(&sum_kk,1,MPI_DOUBLE,0,comm);
  // printf("MPI_Bcast kk myproc=%d, sum_a=%f, sum_kk=%f\n",myproc,sum_a,sum_kk);


// now compute new alphaw for proc=0, and only if variable
if(myproc==0 && prop->alphaw==0){

  // only update if at least 3 good points   
  // alpha(i+1) = alpha(i)+dalpha; 
  if(sum_kk>2){
    average->alphaw+=0.5*(1-sum_a/sum_kk);
  }
  // overwrite any NAN results
  if(IsNan(average->alphaw)){
    average->alphaw=prop->alphaw;
  }
  // set range 0.1 < alphaw < 10 for stability
  if(average->alphaw<0.1){
    average->alphaw=0.1;
  }
  if(average->alphaw>10){
    average->alphaw=10;
  }

  printf("end avgerage->alphaw=%f, dalpha=%f\n",average->alphaw,0.5*(1-sum_a/sum_kk));
}

// send new alphaw to all processors
MPI_Bcast(&average->alphaw,1,MPI_DOUBLE,0,comm);
// printf("end average->alphaw=%f\n", average->alphaw);

}
*/

/*
 * Function: InitSponge
 * Usage: InitSponge(grid,myproc);
 * -------------------------------
 * Apply a sponge layer to all type 2 boundaries.
 *
 */
void InitSponge(gridT *grid, int myproc) {
  int Nb, p1, p2, mark, g1, g2;
  int j, n, NeAll, NpAll;
  REAL *xb, *yb, *xp, *yp, r2;
  char str[BUFFERLENGTH];
  FILE *ifile;

  NeAll = MPI_GetSize(EDGEFILE,"InitSponge",myproc);
  NpAll = MPI_GetSize(POINTSFILE,"InitSponge",myproc);

  xp = (REAL *)SunMalloc(NpAll*sizeof(REAL),"InitSponge");
  yp = (REAL *)SunMalloc(NpAll*sizeof(REAL),"InitSponge");
  rSponge = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"InitSponge");
  rn1 = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"InitSponge");
  rn2 = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"InitSponge");

  REAL xmax, xmin;




  // Read in points on entire grid
  ifile = MPI_FOpen(POINTSFILE,"r","InitSponge",myproc);
  for(j=0;j<NpAll;j++) {
    xp[j]=getfield(ifile,str);
    yp[j]=getfield(ifile,str);
    getfield(ifile,str);
  }
  fclose(ifile);

  // Count number of nonzero boundary markers on entire grid
  ifile = MPI_FOpen(EDGEFILE,"r","InitSponge",myproc);
  Nb = 0;
  for(j=0;j<NeAll;j++) {
    fscanf(ifile, "%d %d %d %d %d",&p1,&p2,&mark,&g1,&g2);
    if(mark==2 || mark==3)
      Nb++;
  }
  fclose(ifile);

  xb = (REAL *)SunMalloc(Nb*sizeof(REAL),"InitSponge");
  yb = (REAL *)SunMalloc(Nb*sizeof(REAL),"InitSponge");

  n=0;
  ifile = MPI_FOpen(EDGEFILE,"r","InitSponge",myproc);
  for(j=0;j<NeAll;j++) {
    fscanf(ifile, "%d %d %d %d %d",&p1,&p2,&mark,&g1,&g2);
    if(mark==2 || mark==3) {
      xb[n]=0.5*(xp[p1]+xp[p2]);
      yb[n]=0.5*(yp[p1]+yp[p2]);
      n++;
    }
  }
  fclose(ifile);  

  // Now compute the minimum distance between the edge on the
  // local processor and the boundary and place this in rSponge.
  xmax=0;
  xmin=INFTY;
  for(j=0;j<grid->Ne;j++) {
    rSponge[j]=INFTY;

    if(grid->xe[j]>xmax){
      xmax=grid->xe[j];
    }
    if(grid->xe[j]<xmin){
      xmin=grid->xe[j];
    }

    for(n=0;n<Nb;n++) {
      r2=pow(xb[n]-grid->xe[j],2)+pow(yb[n]-grid->ye[j],2);
      if(r2<rSponge[j]){
        rSponge[j]=r2;
        rn1[j]=(xb[n]-grid->xe[j])/sqrt(r2);
        rn2[j]=(yb[n]-grid->ye[j])/sqrt(r2);
      }

    }
    rSponge[j]=sqrt(rSponge[j]);
    //    printf("Processor %d: rSponge[%d]=%f\n",myproc,j,rSponge[j]);
  }
  //printf("initsponge xmin=%f, xmax=%f\n", xmin, xmax);
}
  
