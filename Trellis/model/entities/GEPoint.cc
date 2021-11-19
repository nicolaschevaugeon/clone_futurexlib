/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#include "GEPoint.h"
#include "GEdge.h"

int GEPoint::operator == (const GEPoint &pt)
{ return Gent==pt.Gent && Par==pt.Par; }

void GEPoint::checkPar() const
{
  if(Gent){
    AOMD::Range<double> pb = Gent->parBounds(0);
    if(Par < pb.low() || Par > pb.high())
      std::cerr << "GEPoint - par outside range of edge\n";
  }
}

// rather have this function inline, but can't with hp's current compiler
// since then GEPoint.h needs to include GEdge.h, but GEdge.h needs (for
// hp) to include GEPoint.h to instantiate a SSList<GEPoint>
const AOMD::GEntity * GEPoint::gEnt() const
{ return Gent; }

GPoint * GEPoint::clone() const
{ return new GEPoint(*this); }
