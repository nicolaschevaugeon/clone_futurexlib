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
/* 
 * $Id: GEdge.h,v 1.15 2005/02/02 00:08:54 acbauer Exp $
 * Description: 

*/
#ifndef H_GEdge
#define H_GEdge

#include "GEntity.h"

// needed so that hp compiler can do SSList<GEPoint>
#include "GEPoint.h"  
#include "SBlock.h"

template<class T> class SSList;
class GPoint;
class SVector3;
class SPoint2;
class GEdgeUsePair;
class GEdgeUse;

/** A model edge. */
class GEdge : public AOMD::GEntity {
public:
  GEdge(SGModel *model, int tag, GVertex *v0, GVertex *v1);
  ~GEdge() override ;

  TopoType::Value type() const override;
  int dim() const override;

  SSList<GRegion*> regions() const override;
  SSList<GFace*> faces() const override;
  SSList<GEdge*> edges() const override;
  SSList<GVertex*> vertices() const override;

  Logical::Value periodic(int dim=0) const override;
  Logical::Value continuous(int dim=0) const override;

  /** True if the edge is a seam for the given face. */
  virtual int isSeam(GFace *face) const;

  /** Get period of edge if periodic. */
  virtual double period() const;

  /** Get number of edge uses for this edge. */
  int numUses() const;

  /** Get list of edge uses for this edge. */
  SSList<GEdgeUsePair*> uses() const;

  /** Get vertex on given end (0,1) of edge. */
  GVertex* vertex( int n) const;

  int inClosure(GEntity *ent) const override; // if ent is in closure of this edge

  // Geometric ops
  /** Get parameter on edge for given point. */
  virtual double param(const GPoint &pt);

  /** Get the parameter location for a point in space on the edge. */
  virtual double parFromPoint(const AOMD::SPoint3 &) const = 0;

  /** Get the point for the given parameter location. */
  virtual GEPoint point(double p) const = 0;  

  /** Get the closest point on the edge to the given point. */
  virtual GEPoint closestPoint(const AOMD::SPoint3 & queryPoint);

  /** True if the edge contains the given parameter. */
  virtual int containsParam(double pt) const = 0;

  int classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const override;
  virtual int classifyParam(double par, int chkb, GEntity **bent) const;

  /** Get first derivative of edge at the given parameter. */
  virtual SVector3 firstDer(double par) const = 0;

  /** Get nth derivative at the given paramater. */
  virtual void nthDerivative(double param, int n, double *array) const=0;

  /** reparmaterize the point onto the given face. */
  virtual SPoint2 reparamOnFace(GFace *face, double epar,int dir) const = 0;			  
  // these should not be called by anything other than GEnitity dervied
  // classes
  void setVertex(GVertex *v, int which);
  GEdgeUse * addLoop(GLoopUse *l, int dir);
  void removeUse(GEdgeUsePair *);
  void mergeUses(GShell *su);

  void write(std::ostream &out);
  void read(std::istream &in);
  GEdgeUsePair * use(int id);
  // this should be protected, but Sun CC complains about access from derived
  // classes

  virtual GVertex * split(double par)=0;
  virtual SSList<GEPoint> intersect(int fAxis, double fPar, GFace *f) =0;

  // debugging
  void info() const;
protected:
  virtual int getOrderedFaces(SBlock<GFace*> *fs, SBlock<int> *dirs) const;

  void replace(GEdge *enew);

  void splitAtVertex(GVertex *v, GEdge *es, GEdge *ee);

  GVertex *Vertices[2];
  SSList<GEdgeUsePair *> Uses;

};


#endif
