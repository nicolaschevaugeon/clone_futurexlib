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
#ifndef MMIGRATEUTIL_H
#define MMIGRATEUTIL_H

#include "mMesh.h"
#include "mEntity.h"
#include "mAOMD.h"
#include "pmGraphs.h"
#include "mAttachableDataContainer.h"
#include "AOMD_LoadBalancer.h"
#include "mException.h"
#include <cassert>
#include <vector>
#include <set>
#include <list>
#include "ParUtil.h"

namespace AOMD {
bool isEntityOnCurrentProcessor (mEntity *e, int currentProcessor );
void printEntity (mEntity *e);
void CreateMigrationVectors (mMesh *theMesh,
				    int *partitionVector, 
				    int delta_id,
				    int from,
				    int levelMin,
				    std::list<mEntity *> &toDelete, 
				    int nbProcs,
				    int strategy = 0);
void CreateMigrationVectors_oneLevel(mMesh *theMesh,
                                    int *partitionVector,
                                    int delta_id,
                                    int from,
                                    int levelMin,
                                    std::list<mEntity *> &toDelete,
                                    int nbProcs);

void removeEntities(mMesh*, int, std::list<mEntity*>&);
//void getSortedPids_poor_to_rich(mMesh* mesh,std::vector<int>&);
int getPid_withBalance(mMesh* mesh, mEntity* ent, std::vector<int>&);



static void attachIntVector (mEntity *e, bool mark, unsigned int mk,int part)
{
  if(!mark)
    {
      if(e->getData(mk))e->deleteData(mk);
    }
  else
    {
      mAttachableIntVector *avec;
      if(e->getData(mk))
	{
	  avec = (mAttachableIntVector*)e->getData(mk);
	}
      else
	{
	  avec = new mAttachableIntVector;
	  e->attachData(mk,avec);
	}
      for(int i=0;i<avec->v.size();i++)if(avec->v[i] == part)return;
      avec->v.push_back(part);
    }
}

/*
  e is going to be migrated
  markAllSubs tag entities that are going to be migrated
  e is initially an entity that is to be migrated, we look
  after all sub-entities and tag them also.
*/

static void markAllSubs(mEntity *e, bool mark, int part)

{
  int n = e->getLevel();
  
  unsigned int tagMark = AOMD_Util::Instance()->lookupMeshDataId("_dest");

  attachIntVector (e, mark, tagMark, part);

  std::set<mEntity*> theWholeFamily;
  e->getAllEntitiesOfAllDimensions(theWholeFamily);
  for(std::set<mEntity*>::iterator it2 = theWholeFamily.begin(); 
      it2 != theWholeFamily.end();++it2)
    {
      mEntity *sub = *it2;
      attachIntVector (sub, mark, tagMark, part);
    }
}

static void markAllSubs_oneLevel(mEntity *e, bool mark, int part)
{

  int n = e->getLevel();
 
  unsigned int tagMark = AOMD_Util::Instance()->lookupMeshDataId("_dest");
 
  attachIntVector (e, mark, tagMark, part);
 
  std::set<mEntity*> theWholeFamily;
    e->getAllEntitiesOfAllDimensions_oneLevel(theWholeFamily);
  for(std::set<mEntity*>::iterator it2 = theWholeFamily.begin();
      it2 != theWholeFamily.end();++it2)
    {
      mEntity *sub = *it2;
      attachIntVector (sub, mark, tagMark, part);
    }

}

static int getNbDests(mEntity *e)
{
  unsigned int tagMark = AOMD_Util::Instance()->lookupMeshDataId("_dest");
  if(e->getData(tagMark))
    {
      mAttachableIntVector *avec;
      avec = (mAttachableIntVector*)e->getData(tagMark);
      return avec->v.size();
    }
  else
    {
      return 0;
    }
}

static int getDest(mEntity *e, int i)
{
  unsigned int tagMark = AOMD_Util::Instance()->lookupMeshDataId("_dest");
  if(e->getData(tagMark))
    {
      mAttachableIntVector *avec;
      avec = (mAttachableIntVector*)e->getData(tagMark);
      return avec->v[i];
    }
  else
    {
      return -1;
    }
}

static void packVertices (mEntity *e, int *verts, int &k)
{
  for(int i=0;i<e->size(0);i++)verts[k++] = e->get(0,i)->getId();
  //  return;
  if(mEntity *p = e->parent())
    {
      packVertices (p, verts, k);
    } 
}

/*
  I propose here a strategy for load balancing octrees, this strategy has no limitations I think but is NOT
  OPTIMAL. Hypothesis. Under a certain level of refinement, we will not allow coarsening.
  
  Up to this level, all levels of refinement of a cell are on the same process.

  We will do better in the future, this is not so bad anyway.
*/


}  // end of AOMD

#endif
