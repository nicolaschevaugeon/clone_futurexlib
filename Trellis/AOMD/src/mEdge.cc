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
#include "mEdge.h"
#include "mUnits.h"
#include "mEntityContainer.h"
#include "mVertex.h"
#include <cstdio>
#include <cmath>

namespace AOMD {

mEdge::mEdge(mVertex *v1, mVertex *v2, GEntity *classification)
  : mEntity ()
{
  theAdjacencies[0] = new mAdjacencyContainer (2);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theClassification = classification;
  if(!v1 || !v2)return;
  iD = v1->getRAND() + v2->getRAND();
}

mEdge::~mEdge()
= default;

int mEdge::getLevel()const
{
  return 1;
}

int mEdge::getNbTemplates(int what)const
{
  if(what == 0)
    return 2;
  else
    return 0;
}

void mEdge::setVertices(mVertex *v1, mVertex *v2)  
{
  theAdjacencies[0]->clear();
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);
  iD = v1->getRAND() + v2->getRAND();
}

} // end of namespace
