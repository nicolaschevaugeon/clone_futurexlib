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
#include "mVertex.h"
#include "mException.h"
#include "mVector.h"
#include "mIdGenerator.h"
#include <cstdio>

using std::ostream;

namespace AOMD {

mVertex::mVertex(mIdGenerator &theIdGenerator, const Trellis_Util::mPoint &pt, GEntity *classif)
  :mEntity(),p(pt)
{
  iD = theIdGenerator.generateId();
  theClassification = classif;
  RAND = 100000002;
}


mVertex::mVertex(int theId, const Trellis_Util::mPoint &pt, GEntity *classif)
  :mEntity(),p(pt)
{
  iD = theId;
  theClassification = classif;
  RAND = 100000002;
}

int mVertex::getRAND()
{
  if (RAND == 100000002)
    {
      srand(iD);
      RAND = abs(rand()) % 100000000;
    }
  return RAND;
}

void mVertex::deleteId(mIdGenerator &theIdGenerator)
{
  theIdGenerator.addId(iD);
  iD = 0;
}

mVertex::~mVertex()
= default;

int mVertex::getLevel()const
{
  return 0;
}

/*
  compart two vertices using lexicographic order
  EPS is the tolerance
*/

VertexLexicographicLessThan::VertexLexicographicLessThan(double eps)
  : EPS(eps)
{
}

bool VertexLexicographicLessThan :: operator()(mVertex* ent1, mVertex* ent2) const
{
  return ent1->point().lexicographicLessThan (ent2->point(),EPS);
}

} // end of namespace
