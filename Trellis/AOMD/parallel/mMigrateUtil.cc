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
 *  AOMD/parallel/pmMigrateUtil.cc
 *  Created by Eunyoung Seol, on Mon Sep 30 2003, 10:26:37 EDT
 *
 *  File Content: function definition what are used for migration
 *
 *************************************************************************** </i>*/

#include <cstdio>
#include <iostream>
#include <cassert>
#include "mMigrateUtil.h"

#include "AOMD_OwnerManager.h"
#include "AOMD_LoadBalancer.h"
#include "mAttachableDataContainer.h"
#include "mMesh.h"
#include "mEntity.h"
#include "mException.h"
#include "mEdge.h"
#include "mVertex.h"
#include "mRegion.h"
#include "ParUtil.h"
#include "AOMD_Internals.h"
#include "AOMDInternals.h"
#include "AOMD.h"
#include "AOMD_cint.h"

#include "pmGraphs.h"
#include "mAOMD.h"
#include <list>
using std::list;
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

namespace AOMD {

void mEntity::getFamily_oneLevel(std::list<mEntity*> &family)
{
  int n = getLevel();
  // insert self
  family.push_back(this);
  pPList entities;
  mEntity* iadj;
  void *iter = nullptr;
  switch (n)
  {
    case 3: // faces
           entities = R_faces((pRegion)this, 0);
           while (iadj = (mEntity*)PList_next(entities,&iter))
	     family.push_back(iadj);
	   // edges
	   entities = R_edges((pRegion)this, 0);
	   while (iadj = (mEntity*)PList_next(entities,&iter))
	     family.push_back(iadj);
	   // vertices
	   entities = R_vertices((pRegion)this, 0);
	   while (iadj = (mEntity*)PList_next(entities,&iter))
	     family.push_back(iadj);
	     
	   family.unique();
//	   assert (family.size()==15);
	   break;
    case 2:// edges
	   entities = F_edges((pFace)this, 0, (pVertex)nullptr);
	   while (iadj = (mEntity*)PList_next(entities,&iter))
	     family.push_back(iadj);
	   // vertices
	   entities = F_vertices((pFace)this, 1);
	   while (iadj = (mEntity*)PList_next(entities,&iter))
	     family.push_back(iadj);
	     
	   family.unique();
	   //assert (family.size()==15);
    
    case 1:  // faces
//           entities = E_vertices((pEdge)this, 0);
//	   while (iadj = (pEntity)PList_next(entities,&iter))
	   family.push_back(get(0,0));
   	   family.push_back(get(0,1));
	   assert (family.size()==3);
	   break;
    default: break;
  }  
  PList_delete(entities);
  return;
}

void mEntity::getAllEntitiesOfAllDimensions_oneLevel(std::set<mEntity*> &family)
{
  int n = getLevel();
  // insert self
  family.insert(this);
  mEntity* fc;
  mEntity* eg;
  mEntity* vt;
  if (n==3)
  {
    for (int i=0; i<size(2);++i)
    {
      fc=get(2,i);
      family.insert(fc);
      for (int j=0;j<fc->size(1);++j)
      {
        eg = fc->get(1,j);
        if (family.find(eg) == family.end())
           family.insert(eg);
        for (int k=0; k<2; ++k)
        {
          vt = eg->get(0,k);
          if (family.find(vt) == family.end())
           family.insert(vt);
        } // vt
      } // eg
    } // fc
  }
  else if (n==2) //2D mesh
  {
    for (int j=0;j<size(1);++j)
    {
      eg = get(1,j);
      family.insert(eg);
      for (int k=0; k<2; ++k)
      {
         vt = eg->get(0,k);
         if (family.find(vt) == family.end())
           family.insert(vt);
      } // vt
    } // eg
  }
}

void removeEntities(mMesh* mesh, int dim,
                    std::list<mEntity*>& toRemove)
{

}




bool isEntityOnCurrentProcessor (mEntity *e, int currentProcessor )
{
  if (currentProcessor < 0)return true; 
  int N = getNbDests (e);
  for (int i=0;i<N;i++)
    {
      if (currentProcessor == getDest(e,i))return true;
    }
  return false;
}

void printEntity (mEntity *e)
{
  printf ("Processor %d\n",ParUtil::Instance()->rank());
  e->print();
  printf ("%d dests = ",getNbDests(e));
  for (int i=0;i<getNbDests(e);i++)
    {
      printf ("%d ", getDest(e,i));      
    }
  printf("\n");
}
/*static int getDest(mEntity *e, int i)
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
*/

mCommonBdry *mMesh::getInterprocessorCommonBdry( int nbProcs, int *procs , int dim )
{
  int gid = theOwnerManager->find(nbProcs,procs);

  if (nbProcs == 1) return nullptr;
  
  if(gid)
    {      
      for(cbiter it = cbbegin();it != cbend(); ++it)
	{
	  mCommonBdry *g = (*it);
	  if(g->getId() == gid)
	    {
	      return g;
	    }
	}      
    }
  else
    {
      gid = theOwnerManager->create(nbProcs,procs);
      mCommonBdry *newg = getCommonBdry(gid,dim);
      return newg;
    }
  return nullptr;
}

void mMesh::ParallelReclassifyMeshEntities (int n)
{
  // #ifdef PARALLEL  - static partitioning
  //allCommonBdries.clear();
  int parts[1024];
  for(int DIM=0; DIM<=3; DIM++)
    {
      for(iter it = begin(DIM) ; it!= end(DIM) ; ++it)
	{
	  mEntity *e = *it;
	  
	  e->setCommonBdry((mCommonBdry*)nullptr);
	  
	  int nbDests = getNbDests(e);	  
//          cout<<"("<<M_Pid()<<") ParallelReclassify - "<<e->getUid()
//              <<" nbDests = "<<nbDests<<"\n";
	  for(int i=0;i<nbDests;i++)
	    {
	      parts[i] = getDest(e,i);
	    }	      
//          cout<<"("<<M_Pid()<<") ParallelReclassifyMeshEntities "<<e->getUid()
//              <<"(nbDests ="<<nbDests<<endl;
	  e->setCommonBdry(getInterprocessorCommonBdry(nbDests,parts, e->getLevel()));
	}
    }
  //#endif 
  return;
}

// ***********************************************************
void mMesh::setCB_and_collectEntitiesToRemove(int n,
                                 list<mEntity*>& vtToRemove, 
                                 list<mEntity*>& egToRemove,
                                 list<mEntity*>& fcToRemove)
// ***********************************************************
{

  return;
}

void CreateMigrationVectors (mMesh *theMesh,
				    int *partitionVector, 
				    int delta_id,
				    int from,
				    int levelMin,
				    list<mEntity *> &toDelete, 
				    int nbProcs,
				    int strategy)
{
  int *perproc = new int[nbProcs];
  
  int part;
  int k1(0),k2(0);
 
  for(mMesh::iterall it = theMesh->beginall(from) ; it!= theMesh->endall(from) ; ++it)
    {
      mEntity *e = *it;

      k1++;
      for(int i=0;i<nbProcs;i++)perproc[i] = 0;
      if(theMesh->getRefinementLevel(e) == levelMin)
	{
	  int imax = 0;
	  if(strategy == 1 && e->getData(AOMD_Util::Instance()->getId()))
	    {
	      int id = e->getAttachedInt(AOMD_Util::Instance()->getId()) - delta_id;
	      imax = partitionVector[id];
	    }
	  else
	    {
	      list<mEntity *>leaves;
	      e->getLeaves(leaves);
	      for(list<mEntity*>::const_iterator it2 = leaves.begin();it2 != leaves.end();++it2)
		{
		  mEntity *leaf = *it2;
//                  printf("\t(%d) leaf=%s\n",M_Pid(),leaf->getUid());
		  int id = leaf->getAttachedInt(AOMD_Util::Instance()->getId()) - delta_id;
		  part = partitionVector[id];
		  perproc[part]++;
//                  printf("\t(%d) leaf=%s, id=%d, part=%d\n",M_Pid(),leaf->getUid(), id, part);
		}
	      int max = perproc[0];
	      //	      printf("perproc : ");
	      for(int i=1;i<nbProcs;i++)
	      {
	  //		  printf ("%d ",perproc[i]);
	        if(perproc[i]>max)
	        {
	          max = perproc[i];
	          imax = i;
	        }
	      }	    
	      for(list<mEntity*>::const_iterator it2 = leaves.begin();it2 != leaves.end();++it2)
	      {
	        mEntity *leaf = *it2;
	        int id = leaf->getAttachedInt(AOMD_Util::Instance()->getId()) - delta_id;
	        partitionVector[id] = imax;
//                printf("\t(%d) partitionVector[%d]=%d\n", M_Pid(), id, imax);
	      }	
	    }
	  markAllSubs(e,true,imax);
//          cout<<"("<<M_Pid()<<") markAllSubs("<<e->getUid()
//              <<","<<true<<","<<imax<<endl;

	  list<mEntity *>subs;
	  e->getAllSubTree(subs);
	  //  printf("subtree is %d\n",subs.size());
	  for(list<mEntity*>::const_iterator it2 = subs.begin();it2 != subs.end();++it2)
	    {
	      k2++;
	      mEntity *sub = *it2;
	      if(ParUtil::Instance()->rank() != imax)
		toDelete.push_back(sub);
	    }
	}
    }    
  /*
    complete information with a round of communications
  */
  
  if(k1!=k2)
    ParUtil::Instance()->Msg(ParUtil::ERROR," in file %s line %d\n",__FILE__,__LINE__);
  //  ParUtil::Instance()->Msg(ParUtil::WARNING," exchangin starts\n");

  delete [] perproc;
}

void CreateMigrationVectors_oneLevel(mMesh *theMesh,
                                    int *partitionVector,
                                    int delta_id,
                                    int from,
                                    int levelMin,
                                    list<mEntity *> &toDelete,
                                    int nbProcs)
{
  int part; 
  for(mMesh::iterall it = theMesh->beginall(from) ; it!= theMesh->endall(from) ; ++it)
  {      
    mEntity *e = *it;
    int id = e->getAttachedInt(AOMD_Util::Instance()->getId()) - delta_id;
    part = partitionVector[id];
    markAllSubs_oneLevel(e,true,part);
    if (ParUtil::Instance()->rank()!=part)
      toDelete.push_back(e);
  }

}


}  // end of AOMD
