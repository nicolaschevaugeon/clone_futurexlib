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
 *  AOMD/cint/AOMD_Mesh.cc
 *  Created by Seegyoung Seol, on Mon Dec 15 2003, 10:26:37 EDT
 *
 *  File Content: AOMD c-iterface for Mesh
 *
 *************************************************************************** </i>*/
 
#ifndef SIM
#include "AOMD_cint.h"
#include "AOMD_Internals.h"
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
#include "mException.h"
#include "mAOMD.h"
#include "mDebugUtil.h"

#ifdef DMUM
#include "mDmum.h"
#endif

#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cstring>

using namespace AOMD;
using std::cout;
using std::endl;
using std::vector;
using std::copy;
using std::string;

// **********************************
// General mesh Info
// **********************************
int M_getMaxVId(pMesh mesh)
{

  return mesh->getIdGenerator().getMaxValue(); 
 
}

int M_pid()
{
  return P_pid();
}

int M_size()
{
  return P_size();
}

int PM_load (pMesh pm, const char *fName)
{

  return 1;
}
  
int PM_write(pMesh pm, const char *fName)
{

  return 1;
}
  

void M_getCBEntities(pMesh mesh,int dim, std::vector<pEntity>& cbEntities)
{

}

void M_getPOEntities(pMesh mesh, std::vector<pEntity>* poEntities)
{
  for (int i=0; i<3; ++i)
    for (mMesh::iterall it=mesh->beginall(i); it!=mesh->endall(i); ++it)
    {
      if ((*it)->size(i+1)==0)
        poEntities[i].push_back(*it);
    }
  for (mMesh::iterall it=mesh->beginall(3); it!=mesh->endall(3); ++it)
    poEntities[3].push_back(*it);
}

void M_boundingPids(std::vector<int>& bps)
{

}


void M_print(pMesh pm)
{
  ParUtil::Instance()->Msg(ParUtil::INFO,"\n***  M_print(mesh)  ***\n");	  
  pm->printAll();
  ParUtil::Instance()->Barrier(__LINE__,__FILE__);
  ParUtil::Instance()->Msg(ParUtil::INFO,"\n***  M_print(mesh)  ***\n");	  
}


bool M_compare(pMesh mesh, pMesh mesh2)
{
  double t1 = ParUtil::Instance()->wTime();
  ParUtil::Instance()->Msg(ParUtil::INFO,"\n****************************\n");
  ParUtil::Instance()->Msg(ParUtil::INFO,"    M_compare");
  
  mEntity* ent;
  mEntity* ent2;
  pmEntity* pe;
  pmEntity* pe2;
  for (int i=0; i<=3;++i)
  {
  // compare the total number of entities
    assert(mesh->size(i)==mesh2->size(i));
    mMesh::iterall it=mesh->beginall(i);
    mMesh::iterall it2=mesh2->beginall(i);
    mMesh::iterall itend=mesh->endall(i);
    mMesh::iterall itend2=mesh2->endall(i);
    // compare entity by entity
    for (;it!=itend&&it2!=itend2;)
    {
      ent=*it;
      ent2=*it2;
      assert(ent->getId()==ent2->getId());
      assert(ent->getUid().compare(ent2->getUid())==0);
      assert(ent->getLevel()==ent2->getLevel());
      assert(ent->getOwner()==ent2->getOwner());
      assert(ent->getClassification()==ent2->getClassification());

      pe = ent->getPClassification();
      if (pe)
      {
        pe2=ent2->getPClassification();
	assert(pe2!=nullptr);
	assert(pe->getId()==pe2->getId());
	assert(pe->getLevel()==pe2->getLevel());
	assert(pe->getOwner()==pe2->getOwner());
	assert(pe->getNumBPs()==pe2->getNumBPs());
      }
      ++it;
      ++it2;
    }    
  }
  double t2 = ParUtil::Instance()->wTime();
  ParUtil::Instance()->Msg(ParUtil::INFO," SUCCEED\n  ");
  ParUtil::Instance()->Msg(ParUtil::INFO,"   (t = %f sec)\n",t2-t1);
  ParUtil::Instance()->Msg(ParUtil::INFO,"****************************\n\n");
  return true;
}

bool M_verify(pMesh mesh)
{

  double t1 = ParUtil::Instance()->wTime();
  if (ParUtil::Instance()->rank()==0)
  {
    printf("\n****************************\n");
    printf("    M_verify");
  }
  bool ret=verify_mesh(mesh);

  double t2 = ParUtil::Instance()->wTime();
  if (ParUtil::Instance()->rank()==0)
  {
    if(ret) printf(" SUCCEED\n");
    else printf(" FAIL\n");
    printf("   (t = %f sec)\n",t2-t1);
    printf("****************************\n\n");
  }
  return ret;
}

void M_printNumEntities(pMesh mesh)
{
  ParUtil::Instance()->Msg(ParUtil::INFO,"\n****************************\n");
  ParUtil::Instance()->Msg(ParUtil::INFO,"    M_printNumEntities\n\n");
  mesh->printNumEntities();  
  ParUtil::Instance()->Msg(ParUtil::INFO,"\n****************************\n\n");	  
}


int PM_merge(pMesh mesh)
{

  return 0;

}



int M_globalMaxDim(pMesh mesh)
{
  int dim = mesh->getDim();
  int maxDim = P_getMaxInt(dim);
  return maxDim;
}

int M_globalMinDim(pMesh mesh)
{
  int dim = mesh->getDim();
  int maxDim = P_getMinInt(dim);
  return maxDim;
}

void M_numEntitiesOwned(pMesh mesh, int dim, vector<int>& output)
{

  output.push_back(mesh->size(dim));

}

void M_numEntities(pMesh mesh, int dim, vector<int>& output)
{

  output.push_back(mesh->size(dim));

}

void M_updateOwnership(pMesh mesh)
{ 

  return;
}

void DMUM_startMonitoring(pMesh mesh)
{ 
#ifdef DMUM
  mDMUM::Instance()->startMonitor((mMesh*)mesh); 
  if (!P_pid()) system("rm mRep.dat");
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif  
}

void DMUM_stopMonitoring(pMesh mesh)
{
#ifdef DMUM
  mDMUM::Instance()->stopMonitor((mMesh*)mesh); 
  mDMUM::Instance()->write(mesh, "mRep.dat");
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif
}

void DMUM_resumeMonitoring(pMesh mesh)
{  
#ifdef DMUM
  mDMUM::Instance()->resumeMonitor((mMesh*)mesh); 
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif
}

void DMUM_pauseMonitoring(pMesh mesh)
{  
#ifdef DMUM
  mDMUM::Instance()->pauseMonitor((mMesh*)mesh); 
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif
}

void DMUM_print(pMesh mesh)
{
#ifdef DMUM
  mDMUM::Instance()->print((mMesh*)mesh); 
#else 
  if (!P_pid()) cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif
}

bool DMUM_isOn(pMesh mesh)
{
#ifdef DMUM
  return (mMesh*)mesh->DMUM_on; 
#else
  return false;
#endif
}


#endif   /* ifndef SIM */
