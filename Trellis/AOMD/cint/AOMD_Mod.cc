#ifdef SEOLE
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
 *  AOMD/cint/AOMD_Mod.cc
 *  Created by Seegyoung Seol, on Mon Dec 15 2003, 10:26:37 EDT
 *
 *  File Content: AOMD c-iterface for Mesh
 *
 *************************************************************************** </i>*/

#include <iostream>
#include "AOMD.h"
#include "AOMDInternals.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mFace.h"
#include "mRegion.h"
#include "mPoint.h"
#include "mException.h"
#include "mMesh.h"
#include "mAOMD.h"
#include "mBuildAdj.h"

#ifdef DEBUG
using std::cout;
using std::endl;
#endif

using namespace AOMD;

pEdge M_createE(pMesh mesh, pVertex v1, pVertex v2, pGEntity gent)
{
  pEdge p = (pEdge)mesh->createEdge ((mVertex*)v1,(mVertex*)v2,gent);
  v1->add (p);
  v2->add (p);
  return p;
}

pFace M_createF(pMesh mesh, int nEdge, pEdge *edges, int *dirs, pGEntity gent)
{
  pFace f;
  if(nEdge == 3) 
    f = mesh->createFaceWithEdges((mEdge*)edges[0],(mEdge*)edges[1],
 			          (mEdge*)edges[2],gent,dirs);
  else if (nEdge == 4)
    f =  mesh->createFaceWithEdges((mEdge*)edges[0],(mEdge*)edges[1],
			           (mEdge*)edges[2],(mEdge*)edges[3],gent,dirs);
  else 
    return 0;
  for(int i=0;i<nEdge;i++)
    edges[i]->add(f);
  return f;
}

pRegion M_createR(pMesh mesh, int nFace, pFace *faces, int *dirs, pGEntity gent)
{

//  gEntity *ge = convertFromModel (gent,mesh);

  if(nFace == 4) 
    {
      pRegion r = (mRegion*)mesh->createTetWithFaces ((mFace*)faces[0],(mFace*)faces[1],
				(mFace*)faces[2],(mFace*)faces[3],gent);
      faces[0]->add(r);
      faces[1]->add(r);
      faces[2]->add(r);
      faces[3]->add(r);
      return r;
    }
  else return 0;
}

void M_remove(pMesh m, pRegion region, bool do_delete = true) 
{
  int dim = region->getLevel();  
  if(region->isAdjacencyCreated(dim-1))
    {
      for(int i=0;i<region->size(dim-1);i++)
	region->get(dim-1,i)->del(region);  
    }
  if (do_delete)m->DEL(region);
}

void M_removeRegion(pMesh m, pRegion r)
{
  M_remove(m,r);
}

void M_removeFace(pMesh m, pFace face)
{
  M_remove(m,face);
}

void M_removeEdge(pMesh m, pEdge edge)
{
  M_remove(m,edge);
}
#endif
