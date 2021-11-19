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

#ifndef H_TPartialEdge
#define H_TPartialEdge

#include "TEdge.h"

/** An edge that is a portion of another edge. The existing edge
  was split at some given parameter value. */
class TPartialEdge : public TEdge {
public:
  TPartialEdge(SGModel *model, int tag, GEdge *realEdge, GVertex *v0, GVertex *v1, 
	       double parStart, double parEnd);

  AOMD::Range<double> parBounds(int i) const override;
  // Geometric Ops
  GEPoint point(double p) const override;
  GEPoint closestPoint(const AOMD::SPoint3 & queryPoint) override;
  int containsPoint(const AOMD::SPoint3 &pt) const override;
  int containsParam(double par) const override;

  SVector3 firstDer(double par) const override;
  void nthDerivative(double param, int n, double *array) const override;

  SPoint2 reparamOnFace(GFace * face, double epar, int dir) const override;

  GeoRep * geometry() override;

  SSList<GEPoint> intersect(int fAxis, double fPar, GFace *f) override;

protected:
  double parFromPoint(const AOMD::SPoint3 &pt) const override;

  double ParStart, ParEnd;

};

#endif
