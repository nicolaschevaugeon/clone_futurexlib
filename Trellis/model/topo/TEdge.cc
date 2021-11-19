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
#include "TEdge.h"
#include "TPartialEdge.h"
#include "TEdgeSplitVertex.h"
#include <cmath>
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
#include "GRegion.h"
#include "GEPoint.h"
#include "SGModel.h"
#include "SSList.h"
#include "Range.h"

TEdge::TEdge(SGModel *model, int t, GEdge *e,GVertex *v0, GVertex *v1)
: GEdge(model,t, v0,v1), RealEdge(e)
{}

TEdge::~TEdge()
= default;

AOMD::Range<double> TEdge::parBounds(int i) const
{ 
  if (RealEdge)
    return RealEdge->parBounds(i); 
  return AOMD::Range<double>(-1.0e100, 1.0e100);
}

AOMD::SBoundingBox3d TEdge::bounds() const
{ return RealEdge->bounds(); }

GEPoint TEdge::point(double par) const
{ return RealEdge->point(par); }

GEPoint TEdge::closestPoint(const AOMD::SPoint3 & queryPoint)
{ return RealEdge->closestPoint(queryPoint); }

int TEdge::containsPoint(const AOMD::SPoint3 &pt) const
{ return RealEdge->containsPoint(pt); }

int TEdge::containsParam(double par) const
{ return RealEdge->containsParam(par); }

GeoRep * TEdge::geometry()
{ return RealEdge->geometry(); }

SVector3 TEdge::firstDer(double par) const
{ return RealEdge->firstDer(par); }

void TEdge::nthDerivative(double param, int n, double *array) const
{ RealEdge->nthDerivative(param,n,array); }

SPoint2 TEdge::reparamOnFace(GFace * face, double epar, int dir) const
{ return RealEdge->reparamOnFace(face,epar,dir); }

double TEdge::parFromPoint(const AOMD::SPoint3 &pt) const
{ return RealEdge->parFromPoint(pt); }

Logical::Value TEdge::continuous(int i) const
{ return RealEdge->continuous(i); }

Logical::Value TEdge::periodic(int dim) const
{ return RealEdge->periodic(dim); }

Logical::Value TEdge::degenerate(int i) const
{ return RealEdge->degenerate(i); }

int TEdge::isSeam(GFace *face) const
{ return RealEdge->isSeam(face); }

GeomType::Value TEdge::geomType() const
{ return RealEdge->geomType(); }

int TEdge::geomDirection() const
{ return RealEdge->geomDirection(); }

double TEdge::tolerance() const
{ 
  if (RealEdge)
    return RealEdge->tolerance(); 
  return 1.0e-6;
}

GVertex * TEdge::split(double par)
{
  SGModel *m=model();
  int t = RealEdge->tag();
  TEdgeSplitVertex *v = new TEdgeSplitVertex(m,t,RealEdge,par);
  TPartialEdge *el = new TPartialEdge(m,t,RealEdge,vertex(0),v,parBounds().low(),par);
  TPartialEdge *eh = new TPartialEdge(m,t,RealEdge,v,vertex(1),par,parBounds().high());
  model()->add(v);
  model()->add(el);
  model()->add(eh);
  GEdge::splitAtVertex(v,el,eh);
  return v;
}

SSList<GEPoint> TEdge::intersect(int fAxis, double fPar, GFace *f)
{
  SSList<GEPoint> ints = RealEdge->intersect(fAxis,fPar,f);
  SSList<GEPoint> tInts;
  GEPoint ept;
  for(SSListIter<GEPoint> eIter(ints); eIter(ept); )
    tInts.append(GEPoint(ept,this,ept.par()));
  return tInts;
}

AOMD::GEntity * TEdge::getBaseRep()
{ return RealEdge->getBaseRep(); }

void * TEdge::getNativePtr() const
{ return RealEdge->getNativePtr(); }

int TEdge::getNativeInt() const
{ return RealEdge->getNativeInt(); }

//  Attribute * TEdge::attribute(const SString &type) const
//  {
//    if(RealEdge)
//      return RealEdge->attribute(type);
//    return  SModelMember::attribute(type); 
//  }

//  SSList<Attribute *> TEdge::attributes(const SString &type) const
//  {
//    if(RealEdge)
//      return RealEdge->attributes(type);
//    return  SModelMember::attributes(type); 
//  }

//  SSList<Attribute *> TEdge::attributes() const
//  {
//    if(RealEdge)
//      return RealEdge->attributes();
//    return  SModelMember::attributes(); 
//  }
