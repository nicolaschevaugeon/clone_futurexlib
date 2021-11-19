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
#ifndef H_NullEdge
#define H_NullEdge

//#include "SList.h"
#include "GEdge.h"
#include "NullEntity.h"

template<class T> class SSList;
class GVertex;
class NullFace;
class NullVertex;

/** A "fake" model edge. */
class NullEdge : public GEdge, public NullEntity {
public:
  NullEdge(NullModel *m, int tag, GVertex *v0, GVertex *v1);
  NullEdge(NullModel *m, int tag);
  ~NullEdge() override;

  RepType::Value repType() const override;

  Logical::Value continuous(int dim) const override;
  //virtual Logical::Value periodic(int dim);
  Logical::Value degenerate(int dim) const override;

  AOMD::Range<double> parBounds(int i) const override;
  AOMD::SBoundingBox3d bounds() const override;
  // Geometric Ops
  GEPoint point(double p) const override;
  GEPoint closestPoint(const AOMD::SPoint3 & queryPoint) override;
  int containsPoint(const AOMD::SPoint3 &pt) const override;
  int containsParam(double pt) const override;

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
};


#endif
