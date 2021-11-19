/**************************************************************************** 

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

*****************************************************************************/
#ifndef _AOMD_C_INTERNALS_H_
#define _AOMD_C_INTERNALS_H_

/* 
   AOMD C interface 
*/
#include <cstdio>
#ifdef SIM
#include "MeshSim.h"
#include "MeshSimInternal.h"
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#endif

#include <vector>
#include <list>
#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

  int M_modify_representation(pMesh, 
			      int from, 
			      int to, 
			      int with, 
			      int add_or_remove);  
  int M_create_ids(pMesh);  
  int M_remove_ids(pMesh);  
  void M_save_topology(char *, pMesh);
  int M_dimension(pMesh);
  enum mType {VERTEX,EDGE,TRI,QUAD,HEX,PRISM,PYRAMID,TET};
  mType   M_GetElementType   (pEntity pE);
  void    M_GetVertices      (pEntity pE, std::vector<pVertex> &);
  void    M_GetBoundary      (pEntity pE, std::vector<pEntity> &);
  void    M_GetNeighbors     (pEntity pE, pEntity neigh[2]);
  int     R_edgeDir(pRegion , int n);
  int     M_GletClassificationTag(pEntity);
  double  M_SizeMin (pEntity);
  pEntity M_RetrieveRepresentation (int *,pMesh);
  int *   M_BuildRepresentation (pEntity);
  void    M_SetRepresentationToOneLevel(pMesh);

  /// parallel related stuff, does nothing in serial ...
  /// these functions are ALL blocking i.e. they synchronize
  /// all processors at their beginning. They should then be
  /// used by ALL processors.
  void    M_LoadBalance (pMesh);  /// load balancing
  void    M_BdryLinkSetup (pMesh);  /// sets up links between partitions
  /// assign a unique id (attached int with a tag) to all meshe entities
  /// in the list. Id's start at initialId.

  /// returns the maximal vertex id for all partitions
  int     M_MaximumVertexId (pMesh);
  /// returns how many remote entities exists for a mesh entity
  int     M_numRemoteCopies (pMesh, pEntity);
  pPList  M_getRemoteCopies (pMesh, pEntity);
  int     M_getRemotePID (void *);
  int     M_OwnerProc (pMesh,pEntity);

  void V_setWhatInType (pMesh pm, pVertex, int);//  categorize a mesh vertex to a model entity type
  void E_setWhatInType (pMesh pm, pEdge, int);  //  categorize a mesh edge to a model entity type
  void F_setWhatInType (pMesh pm, pFace, int);  //  categorize a mesh face to a model entity type

  /// retrieve some entities from model tag, id and entity dimensions
  void M_getEntitiesClassified (int gEntityDim, 
				int gEntityTag, 
				int mEntityDim, 
				std::list<pEntity> & myList); 

  int E_getWhatInType (pEntity e);
  int E_getWhatInTag (pEntity e);

#ifndef SIM

// 
//    Added by Eunyoung Seol - Oct 05, 2002
//    Modified 		     - Jun 10, 2003
// 

// **************************
// mesh loading
// **************************
  void M_setRepresentationFlag(pMesh, int i);

// **************************
// mesh info
// **************************
  int M_getMaxDim(pMesh);
  int M_NumUniqueVertices(pMesh);
  int M_NumUniqueEdges(pMesh);
  int M_NumUniqueFaces(pMesh);
  int M_NumUniqueRegions(pMesh);
  
// **************************
//  common boundary entity 
// **************************
  pEntity EN_getRemoteCopy(pMesh, pEntity, int);
  void EN_getRemoteCopies(pMesh,pEntity,
                         std::vector<std::pair<pEntity,int> >&);
  void EN_addRemoteCopy(pMesh,pEntity,int,pEntity);
  void EN_deleteRemoteCopy(pMesh,pEntity);
  AOMD::pmEntity* EN_getCommonBdry(pEntity);
  void EN_setCommonBdry(pEntity, AOMD::pmEntity*);
  
// **************************
//  entity info  
// **************************
  void EN_setId(pEntity, int);
  int EN_ownerProc(pEntity); 
  bool EN_onPB(pMesh,pEntity ent);      // return true if the entity is on PB
  bool EN_onCB(pEntity); 
  const char* EN_getUid(pEntity ent);	// return unique id
  int EN_getId(pEntity ent);		// return id
  bool EN_canDeleteAttachedData(pEntity,int);

// **************************
// Mesh migration
// **************************
  // load balancer that uses default loadbalancer callback that AOMD defined
  int M_LoadBalance3(pMesh);  
  
  // load balancer that uses load balancer callback that is given by input
  int M_LoadBalance2(pMesh,AOMD::AOMD_LoadBalancerCallbacks &);

  int M_migrateToOneProc(pMesh, int, std::list<pEntity>,
                   AOMD::AOMD_LoadBalancerCallbacks &,
                   int dim, std::vector<pEntity>&, std::vector<pEntity>&);
  void M_setEntityOwnership(pMesh);
  
// **************************
// conforming meshAdapt-specific
// **************************
  void M_bdryLinkSetupWithMeshMod(pMesh,int,int);
  void M_unifyTaggedEdges(pMesh pm, pMeshDataId tag, std::list<pEntity>&);
  // sum up all tagged edges over all processors
  void M_computeUniqueId(pMesh, pMeshDataId tag,int*, int*);

// **************************
// general parallel util
// **************************
  void M_sync();	//synchronization
  int M_Pid();		// return pid
  int M_NumPE();	// return number of processes
  int P_numPE();
  void P_sync();
  // Unify the information max and min through all processors
  void M_unify(double* max, double* min); 
  void M_MergeArray(std::vector<int>&); 
  int M_getMaxNum(int num);

// **************************
// debug functions
// **************************
  bool M_PrintAdjacencyInfo(pMesh,int,int);	
  void M_PrintEntityInfo(pEntity);
  void M_PrintNumEntitiesInfo(pMesh);

// **************************
// miscellaneous - asked by skocak, May 01, 2002
// **************************
  void M_AssignUniqueRange (pMesh pm, 
		           pMeshDataId tag1,
			   pMeshDataId tag2, 
			   std::vector<pEntity> &l1,
			   std::vector<int> &l2, 
			   int initialId);

#endif   /* ifndef SIM */

#ifdef __cplusplus
}
#endif
#endif
