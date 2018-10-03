/*
 * File: sources.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for sources.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _sources_h
#define _sources_h

#include "grid.h"
#include "phys.h"
#include "met.h"
#include "averages.h"
#include "boundaries.h"
 
REAL **v_coriolis;
REAL *rSponge, *rn1, *rn2, xmin, xmax;


void MomentumSource(REAL **usource, gridT *grid, physT *phys, propT *prop, averageT *average, boundT *bound);
void KurtSource(gridT *grid, physT *phys, propT *prop, MPI_Comm comm);
void HeatSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, metT *met, int myproc, MPI_Comm comm);
void SaltSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, metT *met);
void InitSponge(gridT *grid, int myproc);
// ** not used anymore **
//void ComputeAlphaw(gridT *grid, propT *prop, averageT *average, int myproc, int numprocs, MPI_Comm comm);

#endif
