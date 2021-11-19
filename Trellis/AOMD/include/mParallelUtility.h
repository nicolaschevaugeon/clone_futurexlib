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

/*<i> ****************************************************************************
 *
 *  AOMD/include/mParallelUtility.h
 *  Created by Eunyoung Seol, on Mon Mar 04 2002, 10:26:37 EDT
 *
 *  File Content: definition of the class ParSyncVector
 *                and declaration of some parallel modules
 *
 *************************************************************************** </i>*/

#ifndef MPARALLEL_UTILITY_H
#define MPARALLEL_UTILITY_H

#include "mMesh.h"
#include "mEntity.h"
#include "mEdge.h"
#include "AOMD_OwnerManager.h"
#include <cassert>


using AOMD::mMesh;
using AOMD::mEntity;

namespace AOMD {



/* 
  This function assign the range to each mesh entity.
  IN: vector<mEntity*> entities   contains mesh entities on each processor
      vector<int> numInt          contains entities[i]'s # of integers to assign
      mMesh* theMesh              mesh handle
      unsigned int tag            tag
      int starting_range          starting integral value of the range
  OUT:vector<pair<int,int>> range contains the range of entities[i]
                                  the first contains starting integer
                                  the second contains the # of intergers
  for example, if, in P1, entities=[v1,v2,v3], numInt=[1,2,4]
                   in P2, entities=[v4,v1,v3], numInt=[5,2,4] 
              then,in P1, range=[(1,2),(3,2),(5,4)]
                   in P2, range=[(9,5),(2,2),(5,4)]

output : 
  range[0] is the number of id's for THIS processor (not including the counterparts)
  range[1] is the starting number for this processor
  range[2] is the total number of id's for all processors

*/ 
void assignUniqueRange(mMesh* theMesh,
                       std::vector<mEntity*> & entities,
	               std::vector<int> & numInt,
                       unsigned int startRangeTag,
                       unsigned int endRangeTag,
                       int initialRangeValue); 

void computeUniqueId(mMesh* theMesh,unsigned int tag, int*, int*);
// void bdryLinkSetupWithMeshMod(mMesh* theMesh, int,int);
void unify(double *max,double *min);
void unifyTaggedEdges(mMesh* theMesh, unsigned int tag,
                      std::list<mEntity*> & VEC);

bool canDeleteAttachedData(mEntity*,int);

}  // end of AOMD
#endif 
