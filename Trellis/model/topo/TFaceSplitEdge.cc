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
#include "TFaceSplitEdge.h"
#include "GFace.h"
#include "GEPoint.h"
#include "SVector3.h"
#include "GeoRep.h"
#include "GFPoint.h"
#include "Pair.h"
#include "SPoint2.h"
using std::cerr;

TFaceSplitEdge::TFaceSplitEdge(SGModel *model, int t, GFace *realFace, GVertex *v0, 
			       GVertex *v1, double parStart, double parEnd,
			       int axis, double fParPos, double fParNeg)
  : GEdge(model, t, v0,v1), RealFace(realFace), ParStart(parStart), 
    ParEnd(parEnd), FixedPar(axis)
{ 
  FPar[0] = fParNeg;
  FPar[1] = fParPos;

  if(ParStart >= ParEnd)
    cerr << "TFaceSplitEdge - constructor error";
}

Logical::Value TFaceSplitEdge::continuous(int dim) const
{ return RealFace->continuous(!FixedPar); }

Logical::Value TFaceSplitEdge::periodic(int dim) const
{ return RealFace->periodic(!FixedPar); }

Logical::Value TFaceSplitEdge::degenerate(int dim) const
{ return RealFace->degenerate(!FixedPar); }

AOMD::Range<double> TFaceSplitEdge::parBounds(int i) const
{ return AOMD::Range<double>(ParStart,ParEnd); }

AOMD::SBoundingBox3d TFaceSplitEdge::bounds() const
{ return RealFace->bounds(); }  // a very conservative bounds

GEPoint TFaceSplitEdge::point(double p) const
{ 
  if(p < ParStart || p > ParEnd)
    cerr << "TFaceSplitEdge::point - error\n";
  GFPoint fp = RealFace->point(facePar(p));
  return GEPoint(fp,this,p);
}

GEPoint TFaceSplitEdge::closestPoint(const AOMD::SPoint3 & queryPoint)
{
  cerr << "TFaceSplit::closestPoint - error\n";
  return GEPoint(0,0,0,nullptr,0);
}

int TFaceSplitEdge::containsPoint(const AOMD::SPoint3 &pt) const
{ 
  if(!RealFace->containsPoint(pt))
    return 0;
  // it's on the face, need to see if within tolerance of edge

  // find the face parameters of the point on the seam closest to pt
  SPoint2 fpt = RealFace->parFromPoint(pt);
  // closest point on edge will have same FixedPar parmeter as edge
  fpt[FixedPar] = FPar[0]; 

  AOMD::SPoint3 ept = RealFace->point(fpt);

  if(dist(ept,pt) <= tolerance())
    return 1;

  return 0;
}

int TFaceSplitEdge::containsParam(double par) const
{
  return (par >= ParStart && par <= ParEnd);
}

SVector3 TFaceSplitEdge::firstDer(double par) const
{
  if(par < ParStart || par > ParEnd)
    cerr << "TFaceSplitEdge::firstDer - error\n";
  Pair<SVector3,SVector3> fder = RealFace->firstDer(facePar(par));
  return FixedPar ? fder.second() : fder.first();
}

void TFaceSplitEdge::nthDerivative(double par, int n, double *array) const
{
  if(par < ParStart || par > ParEnd)
    cerr << "TFaceSplitEdge::nthDerivative - error\n";
  cerr << "TFaceSplitEdge::nthDerivative - not implemented\n";
}

SPoint2 TFaceSplitEdge::reparamOnFace(GFace * face, double epar, int dir) const
{
  //if(face!=RealFace)
  //  cerr << "TFaceSplitEdge::reparamOnFace - error\n";
  SPoint2 fpar;
  fpar[FixedPar] = FPar[dir];
  fpar[!FixedPar] = epar;
  return fpar;
}

GeomType::Value TFaceSplitEdge::geomType() const
{ 
  cerr << "TFaceSplit::geomType - error\n";
  return RealFace->geomType(); 
}

GeoRep * TFaceSplitEdge::geometry()
{ 
  cerr << "TFaceSplit::geometry() - error\n";
  return new GeoRep(this); 
}

int TFaceSplitEdge::geomDirection() const
{ 
  //cerr << "TFaceSplit::geomDirection - error\n";
  //return RealFace->geomDirection(); 
  return 1;
}

double TFaceSplitEdge::tolerance() const
{ return RealFace->tolerance(); }

double TFaceSplitEdge::parFromPoint(const AOMD::SPoint3 &pt) const
{ 
  cerr << "TFaceSplit::parFromPoint - error\n";
  SPoint2 p = RealFace->parFromPoint(pt);
  //if(p < ParStart || p < ParEnd)
  //  cerr << "TFaceSplitEdge::parFromPoint - error\n";
  return 0;
}

SPoint2 TFaceSplitEdge::facePar(double p) const
  // this return one of the possibly two face parameters
{
  return SPoint2(FixedPar?p:FPar[1],FixedPar?FPar[1]:p);
}

SSList<GEPoint> TFaceSplitEdge::intersect(int fAxis, double fPar, GFace *f)
{
  notImplemented();
  return SSList<GEPoint>();
}

GVertex * TFaceSplitEdge::split(double par)
{
  notImplemented();
  return nullptr;
}

//  Attribute * TFaceSplitEdge::attribute(const SString &type) const
//  {
//    if(RealFace)
//      return RealFace->attribute(type);
//    return 0;
//  }

//  SSList<Attribute *> TFaceSplitEdge::attributes(const SString &type) const
//  {
//    if(RealFace)
//      return RealFace->attributes(type);
//    return SSList<Attribute*>();
//  }

//  SSList<Attribute *> TFaceSplitEdge::attributes() const
//  {
//    if(RealFace)
//      return RealFace->attributes();
//    return SSList<Attribute*>();
//  }
