/******************************************************************************** 

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written 
   and maintained by the Scientific Computation Research Center (SCOREC) at 
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*********************************************************************************/

/*<i> ****************************************************************************
 *
 *  AOMD/cint/AOMD_cint.h
 *  Created by Seegyoung Seol, on Mon Dec 15 2003, 10:26:37 EDT
 *
 *  File Content: AOMD c-iterface 
 *
 *************************************************************************** </i>*/
 
#ifndef _AOMD_CINT_H_
#define _AOMD_CINT_H_

#include <cstdio>
#include "AOMD.h"
#include "AOMDInternals.h"

#include <vector>
#include <list>
#include <iostream>
#include <map>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef SIM

#ifdef FLEXDB
pPList  AOMD_createAdjPList (pEntity, int);
#endif
/*************************************
  General Parallel operators
**************************************/
double P_getMaxDbl(double);
double P_getMinDbl(double);
double P_getSumDbl(double);
int P_getMaxInt(int num);
int P_getMinInt(int num);
int P_getSumInt(int i);
int P_size();
void P_barrier();
int P_pid();
void P_unifyMinMax(double* max, double* min); 
void P_mergeArray(std::vector<int>&); 
void P_getGlobal(int in, std::vector<int>&);
double P_wtime();
/*************************************
  Mesh operators
**************************************/
int M_getMaxVId(pMesh mesh);
int M_pid();
int M_size();

// the purpose of PM_load/write is not to visualise partitioned mesh,
// but to resume parallel simulation.
int PM_load(pMesh, const char *filename);
int PM_write(pMesh pm, const char *fName);
bool M_compare(pMesh, pMesh);
void M_boundingPids(std::vector<int>& bps);
void M_print(pMesh);
void M_printNumEntities(pMesh);
void M_numEntities(pMesh, int dim, std::vector<int>& );  // asked by Onka
void M_numEntitiesOwned(pMesh, int dim, std::vector<int>& );
void M_updateOwnership(pMesh);
//void M_setEntityOwnership(pMesh);
bool M_verify(pMesh);
int PM_merge(pMesh);

/*************************************
  Entity operators
**************************************/
bool EN_isGhost(pMesh, pEntity);
std::string EN_getUidStr(pEntity);
void EN_print(pEntity);
int EN_owner(pEntity);
bool EN_duplicate(pEntity);
pEntity EN_getCopy(pEntity, int);
void EN_getCopies(pEntity,
              std::vector<std::pair<pEntity,int> >&);
void EN_addCopy(pEntity,int,pEntity);
void EN_clearCopies(pEntity);
AOMD::pmEntity* EN_getPClassification(pEntity);
void EN_setPClassification(pEntity, AOMD::pmEntity*);
void R_center(pEntity, double xyz[3]);
double R_Volume(pEntity);
void F_center(pEntity, double xyz[3]);
void E_center(pEntity, double xyz[3]);
void EN_setWeight(pEntity, double);
int EN_getWeight(pEntity, double*);

/*************************************
  Mesh LB + migration operators
**************************************/




/*************************************
  General Mesh Info operators
**************************************/
int M_globalMaxDim(pMesh);
int M_globalMinDim(pMesh);
// it returns all the mesh entities of given dimension, which are on CB
// (ex)  dim=-1; all entities on CB
//       dim=0: vertices on CB
//       dim=1: edges on CB
//       dim=2: faces on CB
//       dim=3: regions on CB

void M_getCBEntities(pMesh,int dim, std::vector<pEntity>& cbEntities);
void M_getPOEntities(pMesh mesh, std::vector<pEntity>* poEntities);

/*************************************
  parallel meshAdapt specific operators
**************************************/

/*************************************
  parallel Trellis specific operators
**************************************/
void M_assignUniqueRange(pMesh pm, pMeshDataId tag1, pMeshDataId tag2, 
         std::vector<pEntity> &l1, std::vector<int> &l2, int initialId);

/*************************************
  DMUM + FLEXDB
**************************************/
void DMUM_startMonitoring(pMesh);
void DMUM_stopMonitoring(pMesh);
void DMUM_resumeMonitoring(pMesh);
void DMUM_pauseMonitoring(pMesh);
void DMUM_print(pMesh);
bool DMUM_isOn(pMesh);

int M_adjacencyCost(pMesh, int, int);
void M_printAdjacencyCost(pMesh);
void M_printAdjacency(pMesh, int, int);
int M_load_URR(pMesh, const char *, const char*);
int EN_numAdjacency(pEntity, int);
int EN_adjacency(pEntity, int, std::vector<pEntity>&);

#endif  

#ifdef __cplusplus
}
#endif

#endif
