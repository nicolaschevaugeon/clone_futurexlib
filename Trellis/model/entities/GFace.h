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
 /****************************************************************
 * $Id: GFace.h,v 1.14 2005/02/02 00:08:54 acbauer Exp $
 * Description: 

*/
#ifndef H_GFace
#define H_GFace

#include "GEntity.h"
#include "GEPoint.h"
#include "GFaceUse.h"

//#include "DynBitVector.h"
#include "SPoint2.h"
#include "SVector3.h"
#include "Pair.h"

template<class T> class SSList;
class GPoint;
class GEPoint;
class GFPoint;

/** A model face. 
 */
class GFace : public AOMD::GEntity {
public:
  GFace(SGModel *model, int tag, SSList<GEdge*> edges, SSList<int> dirs);
  ~GFace() override;

  TopoType::Value type() const override;
  int dim() const override;

  SSList<GRegion*> regions() const override;
  SSList<GFace*> faces() const override;
  SSList<GEdge*> edges() const override;
  SSList<GVertex*> vertices() const override;
  
  /** Get the face use on the given side of the face. */
  const GFaceUse * use(int dir) const;
  GFaceUse * use(int dir);

  /* Get the region on the given side of the face. */
  GRegion * region(int dir) const;

  int inClosure(GEntity *ent) const override; 

  // geometric ops
  /** Get the location of any parametric degeneracies on the face in the
    given parametric direction. */
  virtual int paramDegeneracies(int dir, double *par) = 0;

  /** Return the point on the face corresponding to the given parameter. */
  virtual GFPoint point(double par1, double par2) const = 0;
  virtual GFPoint point(const SPoint2 &pt) const = 0;

  /** Return the parameter location on the face given a point on the face. */
  SPoint2 param(const GPoint & pt) const;

  /** Return the parmater location on the face given a point in space that
    is on the face. */
  virtual SPoint2 parFromPoint(const AOMD::SPoint3 &) const = 0;

  /** True if the parameter value is interior to the face. */
  virtual int containsParam(const SPoint2 &pt) const = 0;

  /** Period of the face in the given direction. */
  virtual double period(int dir) const;

  int classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const override;
  virtual int classifyParam(const SPoint2 &par, int chkb, GEntity **bent) const;

  /** Return the point on the face closest to the given point. */
  virtual GFPoint closestPoint(const AOMD::SPoint3 & queryPoint);

  /** Return the normal to the face at the given parameter location. */
  virtual SVector3 normal(const SPoint2 &param) const = 0;

  /** Return the first derivate of the face at the parameter location. */
  virtual Pair<SVector3,SVector3> firstDer(const SPoint2 &param) const = 0;
  
  /** Return the nth derivate of the face at the parametric location. */
  virtual double * nthDerivative(const SPoint2 &param, int n, 
				 double *array) const;

  /** Get which side of the face the given use is on. */
  int whichUse(GFaceUse const * use);

  virtual SSList<GEdge*> splitU(double u);
  virtual SSList<GEdge*> splitV(double v);
  virtual SSList<GEdge*> splitSeamU(double ulow, double uhigh);
  virtual SSList<GEdge*> splitSeamV(double vlow, double vhigh);

  // debugging
  void info() const;

  // should only be called internally (by other GEntity classes)
  virtual void insertEdge(GEdge *e);
  virtual void addLoop(SSList<GEdge*> edges, SSList<int> dirs);
  virtual SVector3 inwardSecant(GEdge *e, double ept, int dir);
  
  void write(std::ostream &out);
  void read(std::istream &in);

  void fixPeriodicPar(SPoint2 &pt) const;

  /* true if the surface underlying the face is periodic and we
     need to worry about that. */
  virtual Logical::Value surfPeriodic(int dim) const = 0;

protected:
  GFace(SGModel *model, int tag);  
  // this constructor is for case when face is split by adding an edge
  // these should really be a list of loops
  GFace(SGModel *model, int tag, GLoopUse *frontLoop, GShell *frontShell, 
	GLoopUse *backLoop, GShell *backShell);

  GFaceUse FrontUse;
  GFaceUse BackUse;

};


#endif


