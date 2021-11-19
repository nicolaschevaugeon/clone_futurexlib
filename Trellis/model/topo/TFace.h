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
 /***************************************************************
 * $Id: TFace.h,v 1.11 2005/02/02 00:08:57 acbauer Exp $
 * Description: 
 * Modification 9/15/98 - added the tagId member data that is returned in the
 * case of that RealFace = Null Pointer 

*/

#ifndef H_TFace
#define H_TFace

//#include "SList.h"
#include "GFace.h"

template<class T> class SSList;
class GRegion;

/* A topological face. May represent a copy of an existing face. */
class TFace : public GFace {
public:
  TFace(SGModel *model, int tag, GFace *f,SSList<GEdge*> edges, 
	SSList<int> dirs);
  ~TFace() override;

  GeoRep * geometry() override;

  //Geometric Ops
  AOMD::Range<double> parBounds(int i) const override;
  int paramDegeneracies(int dir, double *par) override;
  AOMD::SBoundingBox3d bounds() const override;
  GFPoint point(double par1, double par2) const override;
  GFPoint point(const SPoint2 &pt) const override;
  GFPoint closestPoint(const AOMD::SPoint3 & queryPoint) override;
  int containsPoint(const AOMD::SPoint3 &pt) const override;
  int containsParam(const SPoint2 &par) const override;

  double period(int dir) const override;

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

  SSList<GEdge*> splitU(double u) override;
  SSList<GEdge*> splitV(double v) override;

  SSList<GEdge*> splitSeamU(double ulow, double uhigh) override;
  SSList<GEdge*> splitSeamV(double vlow, double vhigh) override;

  AOMD::GEntity *getBaseRep() override;
  void * getNativePtr() const override;
  int getNativeInt() const override;

  // these should only be called internalled (by other GEntity classes)
  void insertEdge(GEdge *e) override;

protected:
  Logical::Value surfPeriodic(int dim) const override;

  SSList<GEdge*> splitAtConstParam(int axis, double par1, double par2);

  SPoint2 parFromPoint(const AOMD::SPoint3 &pt) const override;
  void makeNewLoops(GEdge *e, GFaceUse *fu, GFaceUse *fu1, GFaceUse *fu2, GLoopUse **l1, GLoopUse **l2);
    
  GFace *RealFace;
  int d_periodic[2];
};

#endif
