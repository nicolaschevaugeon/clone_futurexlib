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

#ifndef H_GVPoint
#define H_GVPoint

#include "GPoint.h"

class GVertex;

/** A point that is located on a model vertex. */
class GVPoint : public GPoint {
public:
  GVPoint() {}
  GVPoint(double x, double y, double z, const GVertex *v)
    : GPoint(x,y,z) , Gent(v) {}
  GVPoint(const SPoint3 &pt, const GVertex *v)
    : GPoint(pt) , Gent(v) {}
  GVPoint(const GVPoint &pt)
    = default;

  const AOMD::GEntity *gEnt() const override;

  GPoint *clone() const override;
private:  
  const GVertex *Gent;
};


#endif
