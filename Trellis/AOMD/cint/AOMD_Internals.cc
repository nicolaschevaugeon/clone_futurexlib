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
#include "AOMDfwd.h"
#include "AOMD_Defs.h"
#include "ParUtil.h"
#include "mParallelUtility.h"

#include "AOMD_Internals.h"
#include "mVector.h"
#include "mPoint.h"
#include "mVertex.h"

#include "AOMDfwd.h"
#include "AOMD_Defs.h"
#include "mEntity.h"
#include "mMesh.h"
#include "mEdge.h"
#include "modeler.h"


#include <iostream>
#include <list>
#include <vector>

using namespace AOMD;
using std::cout;
using std::endl;
using std::vector;

#ifndef SIM


//*******************************************
int M_dimension(pMesh m)
//*******************************************
{
  return((m->size(3))? 3 : ((m->size(2)) ? 2 : 1 ));
}
#endif
//*******************************************
int M_modify_representation(pMesh pm, int from, int to, int with, int add_or_remove)
//*******************************************
{
#ifndef SIM
  //  printf("changing %d %d %d %d\n",from,to,with,add_or_remove);
  pm->modifyState(from,to,add_or_remove,with);
#endif
  return 1;
}

//*******************************************
void M_SetRepresentationToOneLevel(pMesh pm)
//*******************************************
{
#ifndef SIM     // AOMD related
  if (pm->getRepresentationFlag()==1)
    return;

  double t1=ParUtil::Instance()->wTime();
  ParUtil::Instance()->Msg(ParUtil::INFO,"\n***** M_setRepresentationToOneLevel *****\n");  
  pm->modifyState(3,2,true,0); // region->face
  pm->modifyState(2,1,true,0); // face->edge
  //  pm->modifyState(3,0,false,0);// region->vertex (delete)
  //  pm->modifyState(2,0,false,0);// face->vertex (delete)
  pm->modifyState(2,3,true,0); // face->region
  pm->modifyState(1,2,true,0); // edge->face
  pm->modifyState(0,1,true,0); // vertex -> edge

  double t2=ParUtil::Instance()->wTime();     
  ParUtil::Instance()->Msg(ParUtil::INFO,"\t(t = %f sec)\n", t2-t1);
  ParUtil::Instance()->Msg(ParUtil::INFO,"*****************************************\n");  
  pm->setRepresentationFlag(1);
#else
  throw 1;  
#endif
}

//*******************************************
mType M_GetElementType     (pEntity pE)
//*******************************************
{
#ifndef SIM  	// AOMD related
  return (mType)pE->getType();
#else
  switch(EN_type(pE))
    {
    case 0: 
      return VERTEX;
      break;
    case 1: 
      return EDGE;
      break;
    case 2:
      {
	int numEdges = F_numEdges((pFace)pE);
	if(numEdges == 3)return TRI;
	else if(numEdges == 4)return QUAD;
	else throw 1;	
      }
      break;
    case 3:
      {
	int numFaces = R_numFaces((pRegion)pE);
	if(numFaces == 4)return TET;
	else if(numFaces == 6)return HEX;
	else if(numFaces == 5)return PRISM;
	else throw 1;	
      }
      break;
    default:
      throw 1;
    }
#endif
}

//*******************************************
void  M_GetNeighbors (pEntity pE, pEntity neigh[2])
//*******************************************
{
#ifndef SIM	// AOMD
  int dim = pE->getLevel();
  int n = pE->size(dim+1);
  for(int i=0;i<n;i++)neigh[i] = pE->get(dim+1,i); 
#else
  switch(EN_type(pE))
   { 
    case 0:
      break;
    case 1:      
      break;
    case 2:      
      break;
    default :
      throw 1;
    }
  throw 1;
#endif
}

//*******************************************
void  M_GetVertices (pEntity pE, std::vector<pVertex> &v)
//*******************************************
{
#ifndef SIM	// AOMD
  if(pE->getLevel() == 0)v.push_back(pE);
  else
    {
      int n = pE->size(0);
      for(int i=0;i<n;i++)v.push_back(pE->get(0,i));
    }
#else
  switch(EN_type(pE))
    {
    case 1:
      v.push_back(E_vertex((pEdge)pE,0));
      v.push_back(E_vertex((pEdge)pE,1));
      break;
    case 2:
      { 
	void *iter = 0;
	pEntity iadj; 
	pPList vtxs = F_vertices((pFace)pE,1); 
	while (iadj = (pEntity)PList_next(vtxs,&iter))v.push_back((pVertex)iadj);
	PList_delete(vtxs);
      }
      break;
    case 3:      
      { 
	void *iter = 0;
	pEntity iadj; 
	pPList vtxs = R_vertices((pRegion)pE,1);
	while (iadj = (pEntity)PList_next(vtxs,&iter))v.push_back((pVertex)iadj);
	PList_delete(vtxs);
      }
      break;
    default:
      throw 1;
    }
#endif
}

//*******************************************
void M_GetBoundary      (pEntity pE, std::vector<pEntity> &bnd)
//*******************************************
{
#ifndef SIM	// AOMD
  int dim = pE->getLevel();
  int n = pE->size(dim-1);
  for(int i=0;i<n;i++)
    {
      std::list<pEntity> leaves;
      pE->get(dim-1,i)->getLeaves(leaves);
      //      bnd.insert(bnd.begin(), leaves.begin(), leaves.end());
    }
#else
  throw 1;
#endif
}

#ifdef SIM
static int TEdVt[6][2] = {{0,1},{1,2},{0,2},{0,3},
                         {1,3},{2,3}};
#endif

//*******************************************
int R_edgeDir (pRegion pf, int n)
//*******************************************
{
  int x;
#ifdef SIM
  pPList edgelist = R_edges(pf, 1);
  pEdge edge = (pEdge)PList_item(edgelist, n);
  pVertex v1 = E_vertex(edge, 0);
  pPList vertexlist = R_vertices(pf, 1);
  int indx = TEdVt[n][0];
  pVertex v1_bar = (pVertex)PList_item(vertexlist, indx);
  if(v1 == v1_bar) x = 1;
  else x = -1;
  PList_delete(edgelist);
  PList_delete(vertexlist);
#else 		// AOMD
  mEdge *ed =  (mEdge*)pf->get(1,n);
  mEdge *ed2 = (mEdge*)pf->getTemplate(n,1,0);
  mVertex *v1 = ed->vertex(0);
  mVertex *v1x = ed2->vertex(0);
  if(v1 == v1x)x = 1;
  else x = -1;
  if(ed != ed2) delete ed2; 
#endif
  return x;
}

//*******************************************
double M_SizeMin (pEntity pE)
//*******************************************
{
  Trellis_Util::mPoint p[8];
  std::vector<pVertex> v;
  M_GetVertices (pE, v);

  for(int i=0;i<v.size();i++)
    {
      pPoint pp =  V_point(v[i]); 
      p[i](0) = P_x(pp);
      p[i](1) = P_y(pp);
      p[i](2) = P_z(pp);
    }

  switch (M_GetElementType(pE))
    {
    case VERTEX : return 1.0;    
    case EDGE   : 
      {
	Trellis_Util::mVector v1(p[0],p[1]);
	return sqrt (v1*v1);
      }
    case TRI :
      {
	Trellis_Util::mVector v1(p[0],p[1]);
	Trellis_Util::mVector v2(p[1],p[2]);
	Trellis_Util::mVector v3(p[2],p[0]);
	Trellis_Util::mVector v1v2 = v1 % v2;
	double V = sqrt(v1v2 * v1v2);
	double A1 = sqrt (v1 * v1);
	double A2 = sqrt (v2 * v2);
	double A3 = sqrt (v3 * v3);
	return fabs(V) / (A1+A2+A3);	
      }
    case QUAD :
      {
	Trellis_Util::mVector v1(p[0],p[1]);
	Trellis_Util::mVector v2(p[1],p[2]);
	Trellis_Util::mVector v3(p[2],p[0]);
	Trellis_Util::mVector v1v2 = v1 % v2;
	double V = sqrt(v1v2 * v1v2);
	double A1 = sqrt (v1 * v1);
	double A2 = sqrt (v2 * v2);
	double A3 = sqrt (v3 * v3);
	return fabs(V) / (A1+A2+A3);	
      }
    case TET :
      {
	Trellis_Util::mVector v1(p[1],p[0]);
	Trellis_Util::mVector v2(p[2],p[0]);
	Trellis_Util::mVector v3(p[3],p[0]);
	Trellis_Util::mVector v4(p[1],p[2]);
	Trellis_Util::mVector v5(p[1],p[3]);
	Trellis_Util::mVector v1v2 = v1 % v2;
	Trellis_Util::mVector v1v3 = v1 % v3;
	Trellis_Util::mVector v2v3 = v2 % v3;
	Trellis_Util::mVector v4v5 = v4 % v5;
	double V = v1v2 * v3;
	double A1 = sqrt (v1v2 * v1v2);
	double A2 = sqrt (v1v3 * v1v3);
	double A3 = sqrt (v2v3 * v2v3);
	double A4 = sqrt (v4v5 * v4v5);
	return fabs(V) / (A1+A2+A3+A4);
      }
    default : 
      break;
    }
  return 0.0;
}

pEntity M_RetrieveRepresentation (int *x,pMesh theMesh)
{
#ifdef SIM
  throw 1;
#else
  int nbinfo = x[0];
  int type = x[1];
  switch (type)
    {
    case VERTEX : return (pEntity)theMesh->getVertex (x[2]);
    case EDGE : return (pEntity)theMesh->getEdge (theMesh->getVertex (x[2]),theMesh->getVertex (x[3]));
    case TRI : return (pEntity)theMesh->getTri (theMesh->getVertex (x[2]),theMesh->getVertex (x[3]),
					theMesh->getVertex (x[4]));
    case QUAD : return (pEntity)theMesh->getQuad (theMesh->getVertex (x[2]),theMesh->getVertex (x[3]),
					 theMesh->getVertex (x[4]),theMesh->getVertex (x[5]));
    case TET : return (pEntity)theMesh->getTet  (theMesh->getVertex (x[2]),theMesh->getVertex (x[3]),
					theMesh->getVertex (x[4]),theMesh->getVertex (x[5]));
    case HEX : return (pEntity)theMesh->getHex  (theMesh->getVertex (x[2]),theMesh->getVertex (x[3]),
					theMesh->getVertex (x[4]),theMesh->getVertex (x[5]),
					theMesh->getVertex (x[6]),theMesh->getVertex (x[7]),
					theMesh->getVertex (x[8]),theMesh->getVertex (x[9]));
    default : return nullptr;
    }  
#endif
}

int* M_BuildRepresentation (pEntity pE)
{
  std::vector<pVertex> verts;
  M_GetVertices      (pE, verts);
  int *x = new int [2 + verts.size()];
  x[0] = 2 + verts.size();
  x[1] = M_GetElementType (pE);
  for(int i=0;i<verts.size();i++)x[2+i] = EN_id ((pEntity)verts[i]);
  return x;
}


//// PARALLEL RELATED STUFF

void M_LoadBalance (pMesh pm)
{

}

#ifndef SIM
int M_LoadBalance3(pMesh pm)
{

  return 0;
#endif
}


void M_BdryLinkSetup (pMesh pm)
{
#ifndef SIM
  pm->bdryLinkSetup();
#endif
}

// it returns minimum pid
int M_OwnerProc (pMesh pm, pEntity pe)
{
#ifndef SIM
  return pe->getOwner();
#else
  return 0;
#endif  

}

int M_numRemoteCopies (pMesh pm, pEntity pe)
{
#ifndef SIM	
  return pm->theOwnerManager->count(pe);
#else 		
  return 0;	
#endif		
}

pPList  M_getRemoteCopies (pMesh, pEntity)
{  
  return nullptr;
}

int M_getRemotePID (void *)
{
  return 0;
}

int M_MaximumVertexId (pMesh pm)
{
#ifndef SIM
  int localMaxId = pm->getIdGenerator().getMaxValue();
  int globalMaxId = localMaxId;

  return globalMaxId;

  return 0;
#endif
}

int E_getWhatInTag (pEntity e)
{
#ifndef SIM	
  return GEN_tag(e->getClassification());
#else
  return 0;
#endif
}

int E_getWhatInType (pEntity e)
{
#ifndef SIM	
  return GEN_type(e->getClassification());
#else
  return 0;
#endif
}

void V_setWhatInType (pMesh pm, pVertex e, int x)//  categorize a mesh vertex to a model entity type
{
  throw;
#ifndef SIM	
  e->classify (pm->getGEntity (x, 1));	
#else
  return;
#endif
}

void E_setWhatInType (pMesh pm, pEdge e, int x) //  categorize a mesh edge to a model entity type
{
  throw;
#ifndef SIM
  e->classify(pm->getGEntity (x,1));
#else
  return;
#endif
}

void F_setWhatInType (pMesh pm, pFace e, int x)  //  categorize a mesh face to a model entity type
{
  throw;
#ifndef SIM
  e->classify(pm->getGEntity (x,1));	
  printf ("a face is classified on %d %d %p\n",x,1,e->getClassification());
#else
  return;
#endif
}
