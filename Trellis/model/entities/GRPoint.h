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

#ifndef H_GRPoint
#define H_GRPoint

#include "GPoint.h"

class GRegion;

/** A point in a model region. */
class GRPoint : public GPoint {
  const GRegion *Gent;
public:
  GRPoint() {}
  GRPoint(double x, double y, double z, GRegion *v)
    : GPoint(x,y,z) , Gent(v) {}
  GRPoint(const SPoint3 &pt, GRegion *v)
    : GPoint(pt) , Gent(v) {}
  GRPoint(const GRPoint &pt)
    = default;


  const AOMD::GEntity *gEnt() const override;

  GPoint *clone() const override;

};


#endif
