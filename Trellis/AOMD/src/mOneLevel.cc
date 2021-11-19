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
 *  AOMD/src/mOneLevel.cc
 *  Created by Eunyoung Seol, on Wed Oct 08 2003, 10:26:37 EDT
 *
 *  File Content: function definition what are used for migration
 *
 *************************************************************************** </i>*/
 
#include <iostream>
#include <cassert>

#include "mMesh.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mFace.h"
#include "mTet.h"
#include "mHex.h"
#include "mPrism.h"
#include "mPoint.h"
#include "mException.h"
#include "mBuildAdj.h"
#include "mAOMD.h"

#include <cstdio>

using std::list;
using std::set;
using std::istream;
using std::ostream;
using std::cout;
using std::endl;

namespace AOMD {

void mMesh :: DEL_updateAdj(mEntity *e)
  {
    for (int i=0; i<=getDim();++i)
    {
      if (i==e->getLevel())
        continue;
      if(!e->isAdjacencyCreated(i))
        continue;
      int numAdj = e->size(i);
      for (int j=0; j<numAdj;++j)
        e->get(i,j)->del(e);
    }
    del(e);
    delete e;
  } 



void mMesh::DEL_oneLevel(mEntity* e)
{
  int dim = e->getLevel();
  if(e->isAdjacencyCreated(dim-1))
    {
      for(int i=0;i<e->size(dim-1);i++)
        e->get(dim-1,i)->del(e);
    }
  DEL(e);
}

mEdge* mMesh::createEdge_oneLevel(mVertex *v1, mVertex *v2, GEntity *classif)
{  
  mEdge* p = (mEdge*)createEdge (v1,v2,classif);
  v1->add (p);
  v2->add (p);
  return p;
}

mFace* mMesh::createFaceWithEdges_oneLevel(mEdge *e1 , mEdge *e2, mEdge *e3,
                                      GEntity *classif)
{
  mFace* f;
  f =createFaceWithEdges(e1,e2,e3,classif);
  e1->add(f);
  e2->add(f);
  e3->add(f);
  return f;
}

mTet* mMesh::createTetWithFaces_oneLevel(mFace *f1 , mFace *f2, mFace *f3, mFace *f4,
                                  GEntity *classif)
{
  mTet* r = createTetWithFaces(f1,f2,f3,f4,classif);
  f1->add(r);
  f2->add(r);
  f3->add(r);
  f4->add(r);
  return r;
}

} // end of namespace
