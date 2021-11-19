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
#include "TPartialEdge.h"
#include "GEPoint.h"
#include "SVector3.h"
#include "GeoRep.h"
#include "SPoint2.h"

using std::cerr;

TPartialEdge::TPartialEdge(SGModel *model, int t,GEdge *realEdge, GVertex *v0, 
			   GVertex *v1, double parStart, double parEnd)
  : TEdge(model, t, realEdge, v0,v1), ParStart(parStart), ParEnd(parEnd)
{
  if(ParStart >= ParEnd)
    cerr << "TPartialEdge - parameter error\n";
}


AOMD::Range<double> TPartialEdge::parBounds(int i) const
{ return AOMD::Range<double>(ParStart,ParEnd); }

GEPoint TPartialEdge::point(double p) const
{ 
  if(p < ParStart || p > ParEnd)
    cerr << "TPartialEdge::point - error\n";
  return RealEdge->point(p);
}

GEPoint TPartialEdge::closestPoint(const AOMD::SPoint3 & queryPoint)
{
  GEPoint p = RealEdge->closestPoint(queryPoint);
  if(p.par() < ParStart)
    p = RealEdge->point(ParStart);
  else if(p.par() > ParEnd)
    p = RealEdge->point(ParEnd);
  return p;
}

int TPartialEdge::containsPoint(const AOMD::SPoint3 &pt) const
  // Note: just returns what RealEdge returns
{ return RealEdge->containsPoint(pt); }

int TPartialEdge::containsParam(double par) const
{
  return (par >= ParStart && par <= ParEnd);
}

SVector3 TPartialEdge::firstDer(double par) const
{
  if(par < ParStart || par > ParEnd)
    cerr << "TPartialEdge::firstDer - error\n";
  return RealEdge->firstDer(par);
}

void TPartialEdge::nthDerivative(double par, int n, double *array) const
{
  if(par < ParStart || par > ParEnd)
    cerr << "TPartialEdge::nthDerivative - error\n";
  RealEdge->nthDerivative(par,n,array);
}

SPoint2 TPartialEdge::reparamOnFace(GFace * face, double epar, int dir) const
{
  if(epar < ParStart || epar > ParEnd)
    cerr << "TPartialEdge::nthDerivative - error\n";
  return RealEdge->reparamOnFace(face,epar,dir);
}

GeoRep * TPartialEdge::geometry()
{ return new GeoRep(this); }

double TPartialEdge::parFromPoint(const AOMD::SPoint3 &pt) const
{ 
  double p = RealEdge->parFromPoint(pt);
  if(p < ParStart || p < ParEnd)
    cerr << "TPartialEdge::parFromPoint - error\n";
  return p;
}

SSList<GEPoint> TPartialEdge::intersect(int fAxis, double fPar, GFace *f)
{
  SSList<GEPoint> pts = RealEdge->intersect(fAxis,fPar,f);
  SSList<GEPoint> okpts;
  GEPoint ept;
  for(SSListIter<GEPoint> epIter(pts); epIter(ept); ){
    double epar = ept.par();
    if(epar > ParStart && epar < ParEnd)
      okpts.append(ept);
  }
  return okpts;
}


