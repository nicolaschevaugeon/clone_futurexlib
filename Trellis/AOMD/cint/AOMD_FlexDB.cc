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
 *  AOMD/cint/AOMD_FlexDB.cc
 *  Created by E.S. Seol, on Mon Dec 15 2003, 10:26:37 EDT
 *
 *  File Content: FlexDB interface implementation
 *
 *************************************************************************** </i>*/
#include "AOMD_cint.h"
#include "AOMD_Internals.h"
#include "AOMDfwd.h"
#include "AOMD_Defs.h"
#include "ParUtil.h"
#include "mAOMD.h"
#include "mMesh.h"
#include "AOMD.h"
#include <vector>


#ifdef FLEXDB
#include "mFlexDB.h"
#endif

#ifdef DMUM
#include "mDmum.h"
#endif

#include <iostream>

using namespace AOMD;
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

int M_load_URR(pMesh mesh, const char *filename, const char* MRM_file)
{
#ifdef FLEXDB
#ifdef DMUM
  bool wasOn=false;
  if (mesh->DMUM_on) 
  {
    wasOn=true;
    mesh->DMUM_on=false;
  }
#endif  
  double t1 = ParUtil::Instance()->wTime();
  ParUtil::Instance()->Msg(ParUtil::INFO,"\n********************\n  AOMD with FLEXDB\n");  
  char ext[6];
  strcpy(ext,filename+(strlen(filename)-4));
  if(!strcmp(ext,".sms")) // M_load
  {
    loadSms_URR(mesh,filename,MRM_file); 
    double t2 = ParUtil::Instance()->wTime();
    ParUtil::Instance()->Msg(ParUtil::INFO,"   (t = %f sec)\n", t2-t1);
    ParUtil::Instance()->Msg(ParUtil::INFO,"********************\n");  
  }
  else
  {
    cerr<<"M_load_MRM supports only .sms\n"; throw 1;
  }
#ifdef DMUM
  if (wasOn) mesh->DMUM_on=true;
#endif 
#else // no FLEXDB
  cout<<"Compile AOMD with FLEXDB=1\n";
#endif
  return 1;
}

void M_startMonitoring(pMesh mesh)
{ 
#ifdef DMUM
  mDMUM::Instance()->startMonitor((mMesh*)mesh); 
  if (!P_pid()) system("rm _DMUM_.mrm");
  for (int i=0; i<4; ++i)
    mesh->count_trvs[i] = 0;
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif  
}

void M_stopMonitoring(pMesh mesh)
{
#ifdef DMUM
  mDMUM::Instance()->stopMonitor((mMesh*)mesh); 
  mDMUM::Instance()->write_MRM(mesh, "_DMUM_.mrm");
  cout<<"\n*** Stat. of Mesh Usage ***\n\n";
  cout<<"   # M_VIter "<<mesh->count_trvs[0]<<endl;
  cout<<"   # M_EIter "<<mesh->count_trvs[1]<<endl;
  cout<<"   # M_FIter "<<mesh->count_trvs[2]<<endl;
  cout<<"   # M_RIter "<<mesh->count_trvs[3]<<endl;
  cout<<"\n***************************\n\n";  
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif
}

void M_resumeMonitoring(pMesh mesh)
{  
#ifdef DMUM
  mDMUM::Instance()->resumeMonitor((mMesh*)mesh); 
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif
}

void M_pauseMonitoring(pMesh mesh)
{  
#ifdef DMUM
  mDMUM::Instance()->pauseMonitor((mMesh*)mesh); 
#else 
  if (!P_pid())   cout<<"AOMD ERROR: Compile AOMD with DMUM=1\n";
#endif
}

void M_printMRM(pMesh mesh)
{
#ifdef FLEXDB
  print_MRM((mMesh*)mesh); 
#else
  cout<<"AOMD Info: AOMD is set to the fixed general representation\n";
#endif  
}

bool M_isDMUMon(pMesh mesh)
{
#ifdef DMUM
  return (mMesh*)mesh->DMUM_on; 
#else
  return false;
#endif    
}

int EN_numAdjacency(pEntity ent, int dim)
{
#ifndef FLEXDB
  return ent->size(dim);
#else
  if (ent->getLevel() >= dim || ent->g_mesh->MRM[dim]>0)
    return ent->size(dim);
  vector<mEntity*> adjList;
  EN_adjacency(ent, dim, adjList);
  cout<<"EN_numAdjacency("<<ent->getUid()<<","<<dim
      <<") = "<<adjList.size()<<endl;
  return adjList.size();
#endif
}

#ifdef FLEXDB
void EN_adjacencyRecur(mEntity* ent, int dim, 
                      vector<pEntity>& adjEnts, 
		      mEntity* org_ent, int caller)
{
  pMeshDataId visit_tag = MD_lookupMeshDataId("visited mesh entity");
  int tag_value;
  if (ent->g_mesh->MRM[ent->getLevel()][dim]>0)
  {
    for (int i=0; i<ent->size(dim);++i)
    {
      vector<mEntity*> revAdjEnts;
      if (EN_getDataInt(ent->get(dim,i), visit_tag, &tag_value)) continue;
      EN_adjacency(ent->get(dim,i), org_ent->getLevel(), revAdjEnts);
      if (find(revAdjEnts.begin(),revAdjEnts.end(),org_ent)!=revAdjEnts.end())     
      {
        adjEnts.push_back(ent->get(dim,i));
	EN_attachDataInt(ent->get(dim,i), visit_tag, 1);
      }	
    }	
  }
  // check if local traversal works
  for (int d=0; d<=3; ++d)
  {
    if (d==dim || d==org_ent->getLevel() || d==caller) continue;
    if (ent->g_mesh->MRM[ent->getLevel()][d]>0)
    {
      for (int k=0; k<ent->size(d); ++k)
        EN_adjacencyRecur(ent->get(d,k), dim, adjEnts, org_ent, ent->getLevel());
    }	
  }
}		     
#endif

int EN_adjacency(pEntity ent, int dim, std::vector<pEntity>& adjEnts)
{
#ifndef FLEXDB
  pPList adjList;
  switch (ent->getLevel())
  {
    case 0: switch (dim)
           {
	     case 1: adjList = V_edges(ent); break;
	     case 2: adjList = V_faces(ent); break;
	     case 3: adjList = V_regions(ent); break;
	     default: return 0;
	   }
    case 1: switch (dim)
           {
	     case 0: adjList = E_vertices(ent); break;
	     case 2: adjList = E_faces(ent); break;
	     case 3: adjList = E_regions(ent); break;
	     default: return 0;
	   }
    case 2: switch (dim)
           {
	     case 0: adjList = F_vertices(ent, 1); break;
	     case 1: adjList = F_edges(ent, 1, (pVertex)nullptr); break;
	     case 3: adjList = F_regions(ent); break;
	     default: return 0;
	   }
    case 3: switch (dim)
           {
	     case 0: adjList = R_vertices(ent, 1); break;
	     case 1: adjList = R_edges(ent, 1); break;
	     case 2: adjList = R_faces(ent, 1); break;
	     default: return 0;
	   }
    default: return 0;	   
  }
  mEntity* adjE;
  void* iter;
  while ((adjE = (mEntity*)PList_next(adjList, &iter)))
    adjEnts.push_back(adjE);
  PList_delete(adjList);
  return 1;
#else
  int entLevel = ent->getLevel();
  
  if (ent->g_mesh->MRM[dim][dim]<1)
  {
    cout<<"AOMD WARNING: adj "<<entLevel<<"->"<<dim
        <<"unavailable with the current MRM\n";	
    return 0;	
  }
  if (entLevel >= dim || ent->g_mesh->MRM[entLevel][dim]>0)
  {
    for(int i=0;i<ent->size(dim);i++)
      adjEnts.push_back(ent->get(dim,i));
  }
  else
  {
    for (int d=0; d<=3; ++d)
    {
      if (d==entLevel || d==dim) continue;
      if (ent->g_mesh->MRM[entLevel][d])
      {
        pMeshDataId visit_tag = MD_newMeshDataId("visited mesh entity");
        for (int i=0; i<ent->size(d); ++i)
          EN_adjacencyRecur(ent->get(d,i), dim, adjEnts, ent, entLevel);
        // delete tag attached to adjEnts
	for (vector<mEntity*>::iterator vit=adjEnts.begin();
	    vit!=adjEnts.end(); ++vit)
	  EN_deleteData(*vit, visit_tag);
        // delete tag
	MD_deleteMeshDataId(visit_tag);
        // exit function
	if (adjEnts.size()>0) break;
      }
    }
  } // else 
  if (adjEnts.size()>0) return 1;
  else return 0;
#endif
}

#ifdef FLEXDB
int M_adjacencyCostRecur(pMesh mesh, int i, int j, int k)
{
  if (i==k) return 2;
    
  if (mesh->MRM[i][j]>0) return 1;// stored adj

  // check if local traversal works
  int recur;
  for (int d=0; d<=3; ++d)
  {
    if (i==d || d==j || d==k) continue;
    if (mesh->MRM[i][d]>0)
    {
      recur =  M_adjacencyCostRecur(mesh,d,j,k);
      if (recur < 2) 
        return 1;
    }	
  }
  return 2;
}
#endif

// Return 0: if stored
//  	  1: if local traversal needed
//	  2: if global traversal needed or unavailable
int M_adjacencyCost(pMesh mesh, int i, int j)
{
#ifndef FLEXDB
  // one-level adj
  if (abs(i-j)==1 || i==j) return 0;
  else return 1;
#else
  if (!(mesh->MRM[i][i]>0 && mesh->MRM[j][j]>0)) return 2;

  if (mesh->MRM[i][j]>0) // stored adj
    return 0;

    // check if local traversal works
  int recur;
  for (int d=0; d<=3; ++d)
  {
    if (i==d || d==j) continue;
    if (mesh->MRM[i][d]>0)
    {
      recur =  M_adjacencyCostRecur(mesh,d,j,i);
      if (recur < 2) return 1;
    }	
  }

  return 2;
#endif  
}

void M_printAdjacencyCost(pMesh mesh)
{
  int adjCost[4][4];
  for (int i=0; i<=3; i++)
    for (int j=0; j<=3; ++j)
    {
      adjCost[i][j] = M_adjacencyCost(mesh, i, j);
    }
  cout<<"\n*** Adjacency Cost Table ***\n"
      <<"    (0 stored, 1 local, - global or unavailable)\n";
  for (int i=0; i<=3; i++)
  {
    cout<<"[";
    for (int j=0; j<=3; ++j)
    {
      if (adjCost[i][j]<=1) cout<<adjCost[i][j];
      else cout<<"-";
      if (j==3) cout<<"]";
      else cout<<"  ";
    }
    cout<<"\n";
  }
}

void M_printAdjacency(pMesh mesh, int i, int j)
{
  if (M_adjacencyCost(mesh, i, j)>1)
  {
    cout<<"WARNING: adj "<<i<<"->"<<j<<" unavailable with the current MRM\n";
    return;
  }
    
  cout<<"\n*** PRINT Adjacency From "<<i<<"-dim To "
      <<j<<"-dim ***\n";
  double xyz[3];      
  for (mMesh::iterall it=mesh->beginall(i); it!=mesh->endall(i); ++it)
  {
    if (i==0)
    {
      V_coord(*it, xyz);
      cout<<EN_getUidStr(*it)<<"["<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<"] : ";
    }
    else
      cout<<EN_getUidStr(*it)<<" : ";
    pPList adjList = PList_new();
    switch (i)
    {
      case 0: {
               switch (j)
	       {
	         case 1: adjList = V_edges(*it); break;
		 case 2: adjList = V_faces(*it); break;
		 case 3: adjList = V_regions(*it); break;
	       }
             }
	     break;
      case 1: {
               switch (j)
	       {
	         case 0: PList_append(adjList,(void*)(*it)->get(0,0));
		 	PList_append(adjList,(void*)(*it)->get(0,1));
			break;
		 case 2: adjList = E_faces(*it); break;
		 case 3: adjList = E_regions(*it); break;
	       }
             }
	     break;
      case 2: {
               switch (j)
	       {
	         case 0: adjList = F_vertices(*it,1); break;
		 case 1: adjList = F_edges(*it,1,nullptr); break;
		 case 3: adjList = F_regions(*it); break;
	       }
             }
	     break;      
      case 3: {
               switch (j)
	       {
	         case 0: adjList = R_vertices(*it,1); break;
		 case 1: adjList = R_edges(*it,1); break;
		 case 2: adjList = R_faces(*it,1); break;
	       }
             }
	     break;      
    }
    void* iter=nullptr;
    pEntity ent;
    while ((ent = (mEntity*)PList_next(adjList,&iter)))
      cout<<EN_getUidStr(ent)<<"  ";
    cout<<endl;
    PList_delete(adjList);
   }
   cout<<endl;
}
