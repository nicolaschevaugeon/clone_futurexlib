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
#include <iostream>
#include "mAttachableDataContainer.h"
#include "mMirrorEntity.h"
#include "mVertex.h"
#include "mException.h"
#include "mAOMD.h"

#include <cstdio>

namespace AOMD {

//nico mMirrorVertex::mMirrorVertex(int theId)
//nico   : mVertex (theId,mPoint(0,0,0),0)
//nico {
//nico }

void mMirrorVertex::addCopy (mVertex *e)
{ 
  unsigned int mirrorTag   = AOMD_Util::Instance()->lookupMeshDataId   ("_mirror");
  for(int i=0;i<nbCopies();i++)if(copies[i] == e)return;
  copies.push_back(e);
  mAttachableMirrorVertex *att = new mAttachableMirrorVertex;
  att->e = this;
  e->attachData(mirrorTag,att);
}

void mMirrorVertex::print () const
{
  //nico printf("Mirror Vertex %d \n",getId());
  printf("%d copies\n",nbCopies());
  for(int i=0;i<nbCopies();i++)copies[i]->print();
}

} // end of namespace






