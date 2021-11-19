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

#ifndef H_TFaceSplitEdge
#define H_TFaceSplitEdge

#include "GEdge.h"

class SPoint2;

/** An edge that splits a face at a constant
  parameter value. The edge can lie on a seam in the underlying
  parametric space and thus can have two parameter values w.r.t. the
  face parameter that is fixed. 
  */
class TFaceSplitEdge : public GEdge {
public:
  TFaceSplitEdge(SGModel *model, int tag, GFace *realFace, GVertex *v0, 
		 GVertex *v1, 
		 double parStart, double parEnd, int axis, double fParPos,
		 double fParNeg);


  Logical::Value continuous(int dim) const override;
  Logical::Value periodic(int dim) const override;
  Logical::Value degenerate(int dim) const override;

  AOMD::Range<double> parBounds(int i) const override;
  AOMD::SBoundingBox3d bounds() const override;
  // Geometric Ops
  GEPoint point(double p) const override;
  GEPoint closestPoint(const AOMD::SPoint3 & queryPoint) override;
  int containsPoint(const AOMD::SPoint3 &pt) const override;
  int containsParam(double par) const override;

  SVector3 firstDer(double par) const override;
  void nthDerivative(double param, int n, double *array) const override;

  SPoint2 reparamOnFace(GFace * face, double epar, int dir) const override;

  GeomType::Value geomType() const override;

  GeoRep * geometry() override;
  int geomDirection() const override;

  double tolerance() const override;

  GVertex * split(double par) override;
  SSList<GEPoint> intersect(int fAxis, double fPar, GFace *f) override;

protected:
  double parFromPoint(const AOMD::SPoint3 &pt) const override;
  SPoint2 facePar(double p) const;

  GFace *RealFace;
  double ParStart, ParEnd, FPar[2];
  int FixedPar;  // the parameter value that is fixed along this edge

};

#endif
