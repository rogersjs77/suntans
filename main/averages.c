/*
* Averaging functions
* -------------------
* All functions to calculate and output average quantities 
* 
*/ 
#include "averages.h"

/***********************************************
* Private functions
***********************************************/
/*
 * Function: AllocateAverageVariables()
 * ------------------------------------
 * Allocate memory to the average variable arrays
 *
 */
void AllocateAverageVariables(gridT *grid, averageT **average, propT *prop)
{
  int flag=0, i, j, jptr, ib, Nc=grid->Nc, Ne=grid->Ne, Np=grid->Np, nf, k;


  prop->avgctr=0;
  prop->avgtimectr=0;
  prop->avgfilectr=1*prop->ncfilectr;

  // allocate averageical structure
  *average = (averageT *)SunMalloc(sizeof(averageT),"AllocateAverageVariables");


  (*average)->initialavgfilectr=1*prop->ncfilectr;

  // Allocate 3D arrays
  (*average)->uc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->vc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->w = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->nu_v = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->kappa_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->s = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->rho = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->counter = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");

  // Edge variables
  (*average)->U_F = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->s_F = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->T_F = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");

  // edge velocities
  (*average)->u = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  (*average)->u_avg = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  // (*average)->u_store = (REAL ***)SunMalloc(Ne*sizeof(REAL **),"AllocateAverageVariables");

  // edge velocity variance
  // (*average)->uw_var = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  // (*average)->uw_avg = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  // (*average)->uw_var_avg = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocateAverageVariables");
  // (*average)->uw_var_store = (REAL ***)SunMalloc(Ne*sizeof(REAL **),"AllocateAverageVariables");

  if(prop->calcage)
      (*average)->agec = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");
      (*average)->agealpha = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateAverageVariables");

  // cell-centered averageical variables in plan (no vertical direction)
  (*average)->h = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
  (*average)->h_avg = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
  (*average)->s_dz = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
  (*average)->T_dz = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");

  if(prop->metmodel>0){
    (*average)->Uwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Vwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Tair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Pair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->rain = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->RH = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->cloud = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");

    (*average)->Hs = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Hl = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Hlw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->Hsw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->tau_x = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->tau_y = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
    (*average)->EP = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateAverageVariables");
  }
  
  // for each cell allocate memory for the number of layers at that location
  for(i=0;i<Nc;i++) {
      (*average)->uc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->vc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->w[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateAverageVariables");
      (*average)->s[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->rho[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->nu_v[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
     (*average)->kappa_tv[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
     (*average)->counter[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");
      if(prop->calcage)
         (*average)->agec[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");     	
         (*average)->agealpha[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateAverageVariables");     	

  }
  
  // Allocate edge variables for the number of layers at that location
  for(j=0;j<Ne;j++){
      (*average)->U_F[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->s_F[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->T_F[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");

      (*average)->u[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      (*average)->u_avg[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      // (*average)->u_store[j] = (REAL **)SunMalloc(grid->Nke[j]*sizeof(REAL *),"AllocateAverageVariables");

      // (*average)->uw_var[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      // (*average)->uw_avg[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      // (*average)->uw_var_avg[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocateAverageVariables");
      // (*average)->uw_var_store[j] = (REAL **)SunMalloc(grid->Nke[j]*sizeof(REAL *),"AllocateAverageVariables");
      // now loop over number of stored values
      // for(k=0;k<grid->Nke[j];k++){
      //    // (*average)->u_store[j][k] = (REAL *)SunMalloc(prop->ntaveragestore*sizeof(REAL),"AllocateAverageVariables");
      //    (*average)->uw_var_store[j][k] = (REAL *)SunMalloc(prop->ntaveragestore*sizeof(REAL),"AllocateAverageVariables");
      // }

  }

  // Netcdf write variables
  (*average)->tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"AllocateAverageVariables");
  (*average)->tmpvarE = (REAL *)SunMalloc(grid->Ne*grid->Nkmax*sizeof(REAL),"AllocateAverageVariables");
  (*average)->tmpvarW = (REAL *)SunMalloc(grid->Nc*(grid->Nkmax+1)*sizeof(REAL),"AllocateAverageVariables");
 
}//End of function

/*
 * Function: ZeroAverageVariables() 
 * ---------------------------------
 * Zero the average arrays
 */
void ZeroAverageVariables(gridT *grid, physT *phys, averageT *average, propT *prop ,MPI_Comm comm, int myproc){

    int i,j,k,kk,Nc=grid->Nc,Ne=grid->Ne;

  for(i=0;i<Nc;i++) {
    average->w[i][grid->Nk[i]]=0;
    for(k=0;k<grid->Nk[i];k++) {
      average->uc[i][k]=0;
      average->vc[i][k]=0;
      average->w[i][k]=0;
      average->s[i][k]=0;
      average->T[i][k]=0;
      average->rho[i][k]=0;
      average->nu_v[i][k]=0;
      average->kappa_tv[i][k]=0;
      average->counter[i][k]=0;

      if(prop->calcage){
      	average->agec[i][k]=0;
      	average->agealpha[i][k]=0;
      }
    }
    // 2D cell-centred variables
    average->h[i]=0;
    average->h_avg[i]=0;
    average->s_dz[i]=0;
    average->T_dz[i]=0;
    if(prop->metmodel>0){
	average->Uwind[i]=0;
	average->Vwind[i]=0;
	average->Tair[i]=0;
	average->Pair[i]=0;
	average->rain[i]=0;
	average->RH[i]=0;
	average->cloud[i]=0;
	average->Hs[i]=0;
	average->Hl[i]=0;
	average->Hlw[i]=0;
	average->Hsw[i]=0;
	average->tau_x[i]=0;
	average->tau_y[i]=0;
	average->EP[i]=0;
    }
  }

  for(j=0;j<Ne;j++){
    for(k=0;k<grid->Nke[j];k++){
      average->U_F[j][k]=0;
	    average->s_F[j][k]=0;
	    average->T_F[j][k]=0;
      average->u[j][k]=0;
      // average->uw_var[j][k]=0;
      // average->uw_avg[j][k]=0;
    }
  }
  // only zero these on the first step
  // if(prop->avgtimectr<1){
  //   for(j=0;j<Ne;j++){
  //     for(k=0;k<grid->Nke[j];k++){
  //       average->u_avg[j][k]=0;
  //     }
  //   }

  //   if(prop->restartAvgNC==1){// restart file
  //     ReadAverageVariables(grid,prop,phys,average,myproc,comm);
  //     prop->restartAvgNC==0; // reset so it only loads once
  //   }
  // }

}//End of function

/*
 * Function: UpdateAverageVariables() 
 * ---------------------------------
 * Updates the average arrays by summing the last time step
 *
 */
void UpdateAverageVariables(gridT *grid, averageT *average, physT *phys, metT *met, propT *prop,MPI_Comm comm, int myproc){

    int i,iptr,j,jptr,k,Nc=grid->Nc,Ne=grid->Ne;
    int nc1, nc2;
    REAL flx,dz, sdz,Tdz,theta=prop->theta;
    const REAL V0 = 1e6;
    const REAL V0inv = 1.0/V0;
    REAL uftemp3, u_bar, u_prime, U_bar, U_prime, u_bar_tilde, u_prime_tilde;
    REAL z, depth_face;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    average->w[i][grid->Nk[i]]+=phys->w[i][grid->Nk[i]];
    //sdz=0;
    //Tdz=0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      // Counter is needed for drying of cells
      average->counter[i][k]+=1;

      average->uc[i][k]+=phys->uc[i][k];
      average->vc[i][k]+=phys->vc[i][k];
      average->w[i][k]+=phys->w[i][k];
      average->s[i][k]+=phys->s[i][k];
      average->T[i][k]+=phys->T[i][k];
      average->rho[i][k]+=phys->rho[i][k];
      average->nu_v[i][k]+=phys->nu_tv[i][k];
      average->kappa_tv[i][k]+=phys->kappa_tv[i][k];

      if(prop->calcage){
      	average->agealpha[i][k]+=age->agealpha[i][k];
	      average->agec[i][k]+=age->agec[i][k];
	/* Calculate mean age on the fly
        if(age->agec[i][k]>1e-10){
	    // Calculate the mean online here
	    average->agemean[i][k]+=age->agealpha[i][k]/age->agec[i][k];
	}else{
	    average->agemean[i][k]+=0;
        }
	*/
      }
      //sdz+=phys->s[i][k]*grid->dzzold[i][k];
      //Tdz+=phys->T[i][k]*grid->dzzold[i][k];
    }
    // 2D cell-centred variables
    average->h_avg[i]+=phys->h[i];
    //Instantaneous values (used for scalar budget)
    average->h[i]=phys->h[i];
    //average->s_dz[i]=sdz; 
    //average->T_dz[i]=Tdz; 

    if(prop->metmodel>0 && prop->metmodel<4){
	average->Uwind[i]+=met->Uwind[i];
	average->Vwind[i]+=met->Vwind[i];
	average->Tair[i]+=met->Tair[i];
	average->Pair[i]+=met->Pair[i];
	average->rain[i]+=met->rain[i];
	average->RH[i]+=met->RH[i];
	average->cloud[i]+=met->cloud[i];
	average->Hs[i]+=met->Hs[i]*V0;
	average->Hl[i]+=met->Hl[i]*V0;
	average->Hlw[i]+=met->Hlw[i]*V0;
	average->Hsw[i]+=met->Hsw[i]*V0;
	average->tau_x[i]+=met->tau_x[i];
	average->tau_y[i]+=met->tau_y[i];
	average->EP[i]+=met->EP[i]*phys->s[i][grid->ctop[i]]*V0; // Surface salt flux 
    }
  }

  // compute velocities
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[4];jptr++) {
      j = grid->edgep[jptr]; 

    /* velocity decomposition */
    //u = u_bar + u_prime (separate in time, low pass, high pass)
    //u = U + u_tilde (separate in space, baroclinic, barotropic)
    //u = (U_bar + u_tilde_bar) + (U_prime + u_tilde_prime); (bt mean, bc mean, bt tide, internal waves)
    //u_bar = average->u_avg;
    //u_prime = phys->u - average->u_avg;
    //U_bar = depth_avg(u_bar);
    //U_prime = depth_avg(u_prime);
    //u_tilde_bar = u_bar - U_bar;
    //u_tilde_prime = u_prime - U_prime;

    /* compute depth avg velocity */ 
    U_bar = U_prime = 0;
    depth_face = 0;
    z=0;
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      z-=grid->dz[k]/2; // z on faces is at midpoint
      // Flux needs to be consistent with continuity equation
      // uftemp3 = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k]); // this is u
      uftemp3=phys->u[j][k];
      u_bar = average->u_avg[j][k];
      u_prime = uftemp3 - u_bar;
      U_bar += grid->dzf[j][k]*u_bar;
      U_prime += grid->dzf[j][k]*u_prime;
      depth_face += grid->dzf[j][k];
      z-=grid->dz[k]/2;
    }
    U_bar/=depth_face;
    U_prime/=depth_face;

    REAL tau = prop->ntaverage*prop->dt;
    REAL dt = prop->dt;


    for(k=grid->etop[j];k<grid->Nke[j];k++){
      flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k])*grid->dzf[j][k]*grid->df[j]; 
      average->U_F[j][k] += flx*V0;

      // compute velocity decomposition (edge normal)
      // u = (U_bar + u_tilde_bar) + (U_prime + u_tilde_prime);
      // uftemp3 = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k]);
      uftemp3=phys->u[j][k];
      u_bar = average->u_avg[j][k];
      u_prime = uftemp3 - u_bar;
      u_bar_tilde = u_bar - U_bar;
      u_prime_tilde = u_prime - U_prime; 

      // this is inline filtering Wolfram JPO
      average->u_avg[j][k]=(dt/tau)*uftemp3 + (1-dt/tau)*u_bar;
      

      // average->uw_var[j][k] += pow(u_prime_tilde,2)*V0;
      // average->uw_avg[j][k] += u_prime_tilde*V0;
      //flx = phys->u[j][k]*grid->dzf[j][k]*grid->df[j]; 
	    // Flux needs to be consistent with continuity equation	    
	    //flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k])*dz*grid->df[j]; 
	    //average->U_F[j][k] += flx; 
      //flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k]); 	    
    }
  }

}//End of function

/*
 * Function: UpdateAverageScalars() 
 * ---------------------------------
 * Updates the average scalars (temp,salt) arrays by summing the last time step
 * **This funciton is separated as this step needs to go before the 
 *   tracers are updated for the current time step**
 *
 */
void UpdateAverageScalars(gridT *grid, averageT *average, physT *phys, metT *met, propT *prop,MPI_Comm comm, int myproc){

    int i,iptr,j,jptr,k,Nc=grid->Nc,Ne=grid->Ne;
    int nc1, nc2;
    REAL flx,sdz,Tdz,theta=prop->theta;
    const REAL V0 = 1e6;
    const REAL V0inv = 1.0/V0;

  // Calculate the depth-integrated scalar concentrations
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    average->w[i][grid->Nk[i]]+=phys->w[i][grid->Nk[i]];
    sdz=0;
    Tdz=0;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      sdz+=phys->s[i][k]*grid->dzz[i][k];
      Tdz+=phys->T[i][k]*grid->dzz[i][k];
    }
    //Instantaneous values (used for scalar budget)
    average->s_dz[i]=sdz; 
    average->T_dz[i]=Tdz; 
  }

  /*
   * Compute salinity and temperature fluxes using TVD scheme
   */
  if(prop->TVD  ){

    //Salt
    if(prop->beta>0){
	// Compute the scalar on the vertical faces (for horiz. advection)
	HorizontalFaceScalars(grid,phys,prop,phys->s,phys->boundary_s,prop->TVDsalt,comm,myproc); 
  	for(jptr=grid->edgedist[0];jptr<grid->edgedist[4];jptr++) {
	  j = grid->edgep[jptr]; 
	  for(k=grid->etop[j];k<grid->Nke[j];k++){
	      //See equation 85 in SUNTANS paper
	      flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k])*grid->dzf[j][k]*grid->df[j]; 
	      //if(phys->u[j][k]>0)
	      if(phys->utmp2[j][k]>0)
		average->s_F[j][k]+=phys->SfHp[j][k] * flx * V0;
	      else
		average->s_F[j][k]+=phys->SfHm[j][k] * flx * V0;
	    }
	  }
     }

     //Temperature
     if(prop->gamma>0){
	HorizontalFaceScalars(grid,phys,prop,phys->T,phys->boundary_T,prop->TVDtemp,comm,myproc); 
  	for(jptr=grid->edgedist[0];jptr<grid->edgedist[4];jptr++) {
	  j = grid->edgep[jptr]; 
	  for(k=grid->etop[j];k<grid->Nke[j];k++){
	      flx = (theta*phys->u[j][k] + (1.0-theta)*phys->utmp2[j][k])*grid->dzf[j][k]*grid->df[j]; 
	      //if(phys->u[j][k]>0)
	      if(phys->utmp2[j][k]>0)
		average->T_F[j][k]+=phys->SfHp[j][k] * flx * V0;
	      else
		average->T_F[j][k]+=phys->SfHm[j][k] * flx * V0;
	    }
	  }
      }

  }else{ // No TVD
  
  }//End flux calculation
 /*
   * Compute salinity and temperature fluxes using central-difference
  */

    //Salt and temp
/*
//      for(j=0;j<Ne;j++){
//	for(k=0;k<grid->Nke[j];k++) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
       j = grid->edgep[jptr]; 
       for(k=grid->etop[j];k<grid->Nke[j];k++){

	  flx = phys->u[j][k]*grid->dzf[j][k]*grid->df[j]; 

	  nc1 = grid->grad[2*j];
	  nc2 = grid->grad[2*j+1];
	  if(nc1==-1) nc1=nc2;
	  if(nc2==-1) nc2=nc1;

	  average->s_F[j][k]+= 0.5*(phys->s[nc1][k]+phys->s[nc2][k])*flx;
	  average->T_F[j][k]+= 0.5*(phys->T[nc1][k]+phys->T[nc2][k])*flx;
	}
      }
*/


}
/*
 * Function: ComputeAverageVariables() 
 * ---------------------------------
 * Computes the average by dividing by the number averaging steps "ntaverage"
 * 
 */
void ComputeAverageVariables(gridT *grid, averageT *average, physT *phys, metT *met, int ntaverage, propT *prop){

    int i,j,k,kk,Nc=grid->Nc,Ne=grid->Ne;
    const REAL nt = 1.0/((REAL)ntaverage);
    REAL nt3d;
    const REAL V0 = 1e-6;

  for(i=0;i<Nc;i++) {
    average->w[i][grid->Nk[i]]+=phys->w[i][grid->Nk[i]];
    for(k=0;k<grid->Nk[i];k++) {
      if(average->counter[i][k]>0){
          nt3d = 1.0/(average->counter[i][k]);
      }else{
      	nt3d= 1.0;
      }
     
      average->uc[i][k] *= nt3d;
      average->vc[i][k] *= nt3d;
      average->w[i][k] *= nt3d;
      average->s[i][k] *= nt3d;
      average->T[i][k] *= nt3d;
      average->rho[i][k] *= nt3d;
      average->nu_v[i][k] *= nt3d;
      average->kappa_tv[i][k] *= nt3d;

      if(prop->calcage){
	    average->agec[i][k] *= nt3d;
	    average->agealpha[i][k] *= nt3d;
      }
    }
    // 2D cell-centred variables
    average->h_avg[i] *= nt;
    //average->s_dz[i] *= nt;
    //average->T_dz[i] *= nt;
    if(prop->metmodel>0){
	average->Uwind[i] *= nt;
	average->Vwind[i] *= nt;
	average->Tair[i] *= nt;
	average->Pair[i] *= nt;
	average->rain[i] *= nt;
	average->RH[i] *= nt;
	average->cloud[i] *= nt;
	average->Hs[i] *= nt*V0;
	average->Hl[i] *= nt*V0;
	average->Hlw[i] *= nt*V0;
	average->Hsw[i] *= nt*V0;
	average->tau_x[i] *= nt;
	average->tau_y[i] *= nt;
	average->EP[i] *= nt*V0;
    }
  }

  for(j=0;j<Ne;j++){
    for(k=0;k<grid->Nke[j];k++){
	    average->U_F[j][k] *= nt *V0;
	    average->s_F[j][k] *= nt *V0;
	    average->T_F[j][k] *= nt *V0;

      // average->uw_avg[j][k] *= nt *V0;
      // average->uw_var[j][k] *= nt *V0;
      // average->uw_var[j][k] -= pow(average->uw_avg[j][k],2);
      // compute edge velocity u from flux U_F
      average->u[j][k]=  average->U_F[j][k]/(grid->dzf[j][k]*grid->df[j]);
    }
  }

  // update u_store, last number is the most recent
  // for(j=0;j<Ne;j++){
  //   for(k=0;k<grid->Nke[j];k++){
  //       if(prop->avgtimectr<1){ // first step only
  //         for(kk=0;kk<prop->ntaveragestore;kk++){
  //         // fill values with u
  //         // average->u_store[j][k][kk] = 0;//average->u[j][k];
  //         average->uw_var_store[j][k][kk] = 0;//average->u[j][k];
  //         }
  //       }
  //       else{ // all other steps
  //         for(kk=1;kk<prop->ntaveragestore;kk++){
  //         // rotate through previous stored values
  //         // average->u_store[j][k][kk-1] = average->u_store[j][k][kk];
  //         average->uw_var_store[j][k][kk-1] = average->uw_var_store[j][k][kk];
  //         }
  //       }
  //     // add in new last value
  //     // average->u_store[j][k][prop->ntaveragestore-1] = average->u[j][k];
  //     average->uw_var_store[j][k][prop->ntaveragestore-1] = average->uw_var[j][k];
  //   }
  // }
        
 // update u_avg
  // REAL uftemp, uftemp2;
  // for(j=0;j<Ne;j++){
  //   for(k=0;k<grid->Nke[j];k++){
  //     uftemp=0;
  //     uftemp2=0;
  //     for(kk=0;kk<prop->ntaveragestore;kk++){
  //       // uftemp += average->u_store[j][k][kk];
  //       uftemp2 += average->uw_var_store[j][k][kk];
  //     }
  //     // average->u_avg[j][k]=uftemp/prop->ntaveragestore;
  //     average->uw_var_avg[j][k]=uftemp2/prop->ntaveragestore;
  //     //average->U_Favg[j][k]=average->U_F[j][k];
  //   }
  // }

  //printf("compute average variables ok");
}//End of function

/*
 * Function: SendRecvAverages()
 * ----------------------------
 *  Communicate the average values amongst processors
 *  Only really necessary for writing
 *
 */
void SendRecvAverages(propT *prop, gridT *grid, averageT *average, MPI_Comm comm, int myproc){
    // Communicate 2D variables
    ISendRecvCellData2D(average->h,grid,myproc,comm);
    ISendRecvCellData2D(average->s_dz,grid,myproc,comm);
    ISendRecvCellData2D(average->T_dz,grid,myproc,comm);
    if(prop->metmodel>0){
	ISendRecvCellData2D(average->Uwind,grid,myproc,comm);
	ISendRecvCellData2D(average->Vwind,grid,myproc,comm);
	ISendRecvCellData2D(average->Tair,grid,myproc,comm);
	ISendRecvCellData2D(average->Pair,grid,myproc,comm);
	ISendRecvCellData2D(average->rain,grid,myproc,comm);
	ISendRecvCellData2D(average->RH,grid,myproc,comm);
	ISendRecvCellData2D(average->cloud,grid,myproc,comm);
	ISendRecvCellData2D(average->Hs,grid,myproc,comm);
	ISendRecvCellData2D(average->Hl,grid,myproc,comm);
	ISendRecvCellData2D(average->Hlw,grid,myproc,comm);
	ISendRecvCellData2D(average->Hsw,grid,myproc,comm);
	ISendRecvCellData2D(average->tau_x,grid,myproc,comm);
	ISendRecvCellData2D(average->tau_y,grid,myproc,comm);
	ISendRecvCellData2D(average->EP,grid,myproc,comm);
    }

    // Communicate the 3D cell data
    ISendRecvCellData3D(average->uc,grid,myproc,comm);
    ISendRecvCellData3D(average->vc,grid,myproc,comm);
    ISendRecvCellData3D(average->s,grid,myproc,comm);
    ISendRecvCellData3D(average->T,grid,myproc,comm);
    ISendRecvCellData3D(average->rho,grid,myproc,comm);
    ISendRecvCellData3D(average->nu_v,grid,myproc,comm);
    ISendRecvCellData3D(average->kappa_tv,grid,myproc,comm);
    if(prop->calcage)
	ISendRecvCellData3D(average->agec,grid,myproc,comm);
	ISendRecvCellData3D(average->agealpha,grid,myproc,comm);

    ISendRecvWData(average->w,grid,myproc,comm);

    // Communicate the 3D edge data
    ISendRecvEdgeData3D(average->U_F,grid,myproc,comm);
    ISendRecvEdgeData3D(average->s_F,grid,myproc,comm);
    ISendRecvEdgeData3D(average->T_F,grid,myproc,comm);
    ISendRecvEdgeData3D(average->u,grid,myproc,comm);
    ISendRecvEdgeData3D(average->u_avg,grid,myproc,comm);
    // ISendRecvEdgeData3D(average->uw_var,grid,myproc,comm);
    // ISendRecvEdgeData3D(average->uw_var_avg,grid,myproc,comm);
     
}//End of function


/*
 * Function: ReadAverageVariables
 * Usage: ReadPhysicalVariables(grid,phys,prop,myproc,comm);
 * ---------------------------------------------------------
 * This function reads in average variables for a restart run
 * from the restart file defined by prop->StartFID.
 *
 */
void ReadAverageVariables(gridT *grid, propT *prop, physT *phys, averageT *average, int myproc, MPI_Comm comm) {

  int i,iptr,j,jptr,k,Nc=grid->Nc,Ne=grid->Ne;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  // open files

  // *** right now this only reads in u_avg, but other variables
  // could be added here later. ***
  
  // intialize u_avg=0
   for(j=0;j<Ne;j++){
      for(k=0;k<grid->Nke[j];k++){
        average->u_avg[j][k]=0;
      }
    }

  // if a restart file, load in data from file
  if(prop->restartAvgNC==1){// restart file
    prop->restartAvgNC==0; // reset so it only loads once

  if(VERBOSE>1 && myproc==0) printf("Reading from rstore AVG ...\n");

  // MPI_GetFile(filename,DATAFILE,"StartFile","OpenFiles",myproc);
  // sprintf(str,"%s.%d",filename,myproc);
  // prop->StartFID = MPI_FOpen(str,"r","OpenFiles",myproc);
  OpenFiles(prop,myproc);

  // do exact same routine as physio.c to unwind start file

  if(fread(&(prop->nstart),sizeof(int),1,prop->StartFID) != 1)
    printf("Error reading prop->nstart\n");
    // Initialise the netcdf time
  // 
  if(fread(phys->h,sizeof(REAL),grid->Nc,prop->StartFID) != grid->Nc)
    printf("Error reading phys->h\n");
  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->Cn_U[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j])
      printf("Error reading phys->Cn_U[j]\n");
  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->Cn_U2[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j]) //AB3
      printf("Error reading phys->Cn_U2[j]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_W[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_W[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_W2[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_W[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_R[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_R[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->Cn_T[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->Cn_T[i]\n");

  if(prop->turbmodel>=1) {
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->Cn_q[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->Cn_q[i]\n");
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->Cn_l[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->Cn_l[i]\n");

    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->qT[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->qT[i]\n");
    for(i=0;i<grid->Nc;i++) 
      if(fread(phys->lT[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
        printf("Error reading phys->lT[i]\n");
  }
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->nu_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->nu_tv[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->kappa_tv[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->kappa_tv[i]\n");

  for(j=0;j<grid->Ne;j++) 
    if(fread(phys->u[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j])
      printf("Error reading phys->u[j]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->w[i],sizeof(REAL),grid->Nk[i]+1,prop->StartFID) != grid->Nk[i]+1)
      printf("Error reading phys->w[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->q[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->q[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->qc[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->qc[i]\n");

  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->s[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->s[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->T[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->T[i]\n");
  for(i=0;i<grid->Nc;i++) 
    if(fread(phys->s0[i],sizeof(REAL),grid->Nk[i],prop->StartFID) != grid->Nk[i])
      printf("Error reading phys->s0[i]\n");

  // printf("read uavg, mpyroc=%d, FID=%d\n",myproc,prop->StartFID);

  for(j=0;j<grid->Ne;j++) 
    if(fread(average->u_avg[j],sizeof(REAL),grid->Nke[j],prop->StartFID) != grid->Nke[j])
      printf("Error reading average->u_avg[j]\n");


  fclose(prop->StartFID);
  // if(VERBOSE>1 && myproc==0) printf("done reading\n");
  }
}
