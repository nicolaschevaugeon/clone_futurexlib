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

/*****************************************************************************
 *
 *  cint/PAOMD_Internals.cc
 *  Created by Eunyoung Seol, on Thur Jun 12 2003, 10:26:37 EDT
 *
 *  C interface definition for PAOMD
 *
 ****************************************************************************/
 
#ifndef SIM

#include "AOMD_Internals.h"
#include "AOMD_cint.h"
#include "AOMDfwd.h"
#include "AOMD_Defs.h"
#include "ParUtil.h"
#include "mParallelUtility.h"
#include "mEntity.h"
#include "mMesh.h"
#include "mEdge.h"
#include "mVector.h"
#include "mPoint.h"
#include "mVertex.h"


#include <iostream>
#include <list>
#include <vector>

using namespace AOMD;
using std::cout;
using std::endl;
using std::vector;

// **********************************
// mesh loading
// **********************************


void M_setRepresentationFlag(mMesh* mesh, int i)
{
  mesh->setRepresentationFlag(i);
}


int M_getMaxDim(pMesh mesh){
  int dim = mesh->getDim();
  int maxDim = M_getMaxNum(dim);  
  return maxDim;
}
  
  int M_NumUniqueVertices(pMesh pmesh)
{
 return pmesh->numUniqueEntities(0);
}

int M_NumUniqueEdges(pMesh pmesh)
{
 return pmesh->numUniqueEntities(1);
}  

int M_NumUniqueFaces(pMesh pmesh)
{
 return pmesh->numUniqueEntities(2);
}  

int M_NumUniqueRegions(pMesh pmesh)
{
 return pmesh->numUniqueEntities(3);
}  

// **********************************
// common boundary entity
// **********************************


pEntity EN_getRemoteCopy(pMesh mesh,pEntity ent,int pid)
{

  return (pEntity)nullptr;

}

void EN_getRemoteCopies(pMesh mesh,pEntity ent,
                        std::vector<std::pair<pEntity,int> >& remoteCopies)
{

}

void EN_addRemoteCopy(pMesh mesh, pEntity ent, int remotePid, pEntity remoteCopy)
{

}
  
void EN_deleteRemoteCopy(pMesh mesh, pEntity ent)
{

}

AOMD::pmEntity* EN_getCommonBdry(pEntity ent)
{
  return EN_getPClassification(ent);
}

void EN_setCommonBdry(pEntity ent, AOMD::pmEntity* cb)
{
//  assert(ent->getLevel()!=3);
  EN_setPClassification(ent,cb);
}

// **********************************
// entity info
// **********************************
int EN_ownerProc(pEntity ent)
{
  return ent->getOwner();
}
  
bool EN_onPB(pMesh pmesh,pEntity ent)		// return true if the entity is on PB
{

  return false;

}

bool EN_onCB(pEntity ent)		// return true if the entity is on PB
{
  return EN_duplicate(ent);
  
}

const char* EN_getUid(pEntity ent)	// return unique id
{
  return ent->getUid().data();
}

int EN_getId(pEntity ent)		// return id
{
  return ent->getId();
}


void EN_setId(pEntity ent, int id)
{
  ent->setId(id);
}

bool EN_canDeleteAttachedData(pEntity ent,int proc)
{

 return false;

}
// **********************************
// Mesh migration
// **********************************
void M_setEntityOwnership(pMesh pm)
{

  return;
}

int M_LoadBalance2(pMesh pm,AOMD_LoadBalancerCallbacks& cb)
{

  return 0;

}


int M_migrateToOneProc(pMesh pm,int dimToMove,std::list<mEntity*> entList,
                    AOMD_LoadBalancerCallbacks &cb,
                    int dimToKnow,
                    std::vector<pEntity>&rmE,
                    std::vector<pEntity>& newE)
{

  return 0;

}

// **********************************
// conforming meshAdapt-specific
// **********************************
void M_bdryLinkSetupWithMeshMod(pMesh pm, int min, int max)
{

}


void M_unifyTaggedEdges(pMesh mesh, pMeshDataId tag, std::list<pEntity>& edges)
{

}


void M_computeUniqueId(pMesh mesh, pMeshDataId tag,int* min, int* max)
{

}

// **********************************
//general parallel util
// **********************************
void M_sync()		//synchronization
{
  P_barrier();
}

int M_Pid()
{
  return P_pid();
}

int M_NumPE()
{
  return P_size();
}

void P_sync()           //synchronization
{
  P_barrier();
}
 
int P_numPE()
{
  return P_size();
}
// Unify the information max and min through all processors
void M_unify(double* max, double* min)
{

}

void M_MergeArray(std::vector<int>& vec)
{

} 

int M_getMaxNum(int num)
{
  if (M_NumPE()==1) return num;
  throw;
}

// **********************************
// debug functions
// **********************************

  
  bool M_PrintAdjacencyInfo(pMesh pm,int dim,int adjDim)
  {
   for (mMesh::iterall it = pm->beginall(dim); it!=pm->endall(dim);++it)
   { 
      pEntity e = *it;
      int numAdjEnt = e->size(adjDim);     
      cout<<"("<<M_Pid()<<") M_Test: "<<e->getUid()<<": ";
      for (int i=0;i<numAdjEnt;++i)
      {  
        pEntity adjEnt = e->get(adjDim,i);
        if (!adjEnt) 
        { cout<<"ERROR! ";
          cout<<"("<<M_Pid()<<") "<<i<<"'th adjEnt is not accessible\n";
          return false;
        }
        else
          cout<<adjEnt->getUid()<<", "; 
      } 
      cout<<endl;
   }
   cout<<"("<<M_Pid()<<") M_TestAdjacency("<<dim<<","<<adjDim<<") DONE!\n"; 
   return true;
}

// this will print the mesh entity information 
void M_PrintEntityInfo(pEntity pe)
{
  pe->print();
}

// this will print all the mesh entities information including 
// 	- classification
// 	- common boundary

// Note this funtion cannot be called for test/pAOMD3/main.cc
void M_PrintNumEntitiesInfo(pMesh theMesh)
{
  theMesh->printNumEntities();
}  

// **********************************
// miscellaneous - asked by skocak, May 01, 2002
// **********************************
void M_AssignUniqueRange (pMesh pm, 
	  		       pMeshDataId tag1,
			       pMeshDataId tag2, 
			       std::vector<pEntity> &l1,
			       std::vector<int> &l2, 
			       int initialId)
{

  return;
}

#endif   /* ifndef SIM */
