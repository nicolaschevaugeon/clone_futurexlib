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
#ifndef H_NullFace
#define H_NullFace

 /***************************************************************
 * $Id: NullFace.h,v 1.10 2005/02/02 00:08:56 acbauer Exp $
 * Description: 

*/
//#include "SList.h"
#include "GFace.h"
#include "NullEntity.h"

template<class T> class SSList;
class NullEdge;
class NullRegion;
class GRegion;
class GEPoint;

/** A "fake" model face. 
 */
class NullFace : public GFace, public NullEntity {
public:
  NullFace(NullModel *m, int tag, const SSList<GEdge*> &edges, const SSList<int> &dirs);
  NullFace(NullModel *m, int tag);
  ~NullFace() override;

  RepType::Value repType() const override;

  GeoRep * geometry() override;

  //Geometric Ops
  AOMD::Range<double> parBounds(int i) const override;
  int paramDegeneracies(int dir, double *par) override;
  AOMD::SBoundingBox3d bounds() const override;
  GFPoint point(double par1, double par2) const override;
  GFPoint point(const SPoint2 &pt) const override;
  GFPoint closestPoint(const AOMD::SPoint3 & queryPoint) override;

  int containsPoint(const AOMD::SPoint3 &pt) const override;
  int containsParam(const SPoint2 &pt) const override;

  SVector3 normal(const SPoint2 &param) const override;
  Pair<SVector3,SVector3> firstDer(const SPoint2 &param) const override;
  double * nthDerivative(const SPoint2 &param, int n, 
				 double *array) const override;

  GeomType::Value geomType() const override;
  int geomDirection() const override;

  Logical::Value continuous(int dim) const override;
  Logical::Value periodic(int dim) const override;
  Logical::Value degenerate(int dim) const override;

  double tolerance() const override;

  SVector3 inwardSecant(GEdge *e, double ept, int dir) override;

protected:
  Logical::Value surfPeriodic(int dim) const override;

  virtual SPoint2 parFromEdgePar(GEdge *edge, double epar, int dir) const;
  SPoint2 parFromPoint(const AOMD::SPoint3 &) const override;

};

#endif
