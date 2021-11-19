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

#ifndef H_GFPoint
#define H_GFPoint

#include "GPoint.h"
#include "SPoint2.h"

class GFace;

/** A point on a model face. Has a parameter location as well
  as an xyz location.
*/
class GFPoint : public GPoint {
public:
  GFPoint() {}
  GFPoint(double x, double y, double z, const GFace *e, const SPoint2 & par)
    : GPoint(x,y,z), Gent(e), Par(par) {}
  GFPoint(const SPoint3 &pt, const GFace *e, const SPoint2 &par)
    : GPoint(pt), Gent(e), Par(par) {}
  GFPoint(const GFPoint &pt)
    = default;

  const AOMD::GEntity *gEnt() const override;

  GPoint *clone() const override;
  /** The parameter location of the point. */
  SPoint2 par() const
    { return Par; }

  void setParam(double r, double s);
private:
  const GFace *Gent;
  SPoint2 Par;

};


#endif
