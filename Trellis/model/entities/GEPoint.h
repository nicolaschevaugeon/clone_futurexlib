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

#ifndef H_GEPoint
#define H_GEPoint

#include "GPoint.h"

class GEdge;

/** A point evaluated on a model edge. */
class GEPoint : public GPoint {
public:
  GEPoint() {}
  GEPoint(double x, double y, double z, const GEdge *e, double par);
  GEPoint(const SPoint3 &pt, const GEdge *e, double par);
  GEPoint(const GEPoint &pt)
    = default;

  const AOMD::GEntity *gEnt() const override;

  GPoint *clone() const override;
  /** Get parametric location on edge of point. */
  double par() const;
  void setPar(double par);

  int operator == (const GEPoint &pt);
private:
  void checkPar() const;

  const GEdge *Gent;
  double Par;

};

inline GEPoint::GEPoint(double x, double y, double z, const GEdge *e, double par)
  : GPoint(x,y,z), Gent(e), Par(par) 
{ 
#ifdef DEBUG
  checkPar(); 
#endif
}

inline GEPoint::GEPoint(const SPoint3 &pt, const GEdge *e, double par)
  : GPoint(pt), Gent(e), Par(par) 
{ 
#ifdef DEBUG
  checkPar(); 
#endif
} 

inline double GEPoint::par() const
{ return Par; }

inline void GEPoint::setPar(double par) 
{ 
  Par = par; 
#ifdef DEBUG
  checkPar(); // debug
#endif
}


#endif
