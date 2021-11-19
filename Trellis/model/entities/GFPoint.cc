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
#include "GFPoint.h"
#include "GFace.h"

void GFPoint::setParam(double r, double s)
{ 
  Par[0] = r;
  Par[1] = s;
  
  SPoint3 loc = Gent->point(r,s);
  P[0] = loc[0];
  P[1] = loc[1];
  P[2] = loc[2];
}

// see note in GEPoint
const AOMD::GEntity * GFPoint::gEnt() const
{ return Gent; } 

GPoint * GFPoint::clone() const
{ return new GFPoint(*this); }
