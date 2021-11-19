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
#include "NullEdge.h"
#include "NullVertex.h"
#include <cmath>
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
//#include "GShell.h"
#include "GRegion.h"
#include "GEPoint.h"
#include "SSList.h"
#include "NullGeoRep.h"
#include "NullModel.h"
#include "NullFace.h"

#include <cfloat>
#define MAXDOUBLE DBL_MAX

NullEdge::NullEdge(NullModel *m, int t, GVertex *v0, GVertex *v1) 
: GEdge(m, t, v0,v1), NullEntity() 
{}

NullEdge::NullEdge(NullModel *m, int t) 
: GEdge(m, t, nullptr,nullptr), NullEntity() 
{}

NullEdge::~NullEdge()
= default;

RepType::Value NullEdge::repType() const
{
  return RepType::Unknown;
}

AOMD::Range<double> NullEdge::parBounds(int i) const
{ 
  notImplemented();
  return( AOMD::Range<double>(0,0) );
}

AOMD::SBoundingBox3d NullEdge::bounds() const
{
  AOMD::SBoundingBox3d box;
  notImplemented();
  return box;
}

GEPoint NullEdge::point(double par) const
{
  notImplemented();
  return GEPoint(AOMD::SPoint3(0,0,0),nullptr,par); // invalid parameter
}

GEPoint NullEdge::closestPoint(const AOMD::SPoint3 & queryPoint)
{
  notImplemented();
  return GEPoint(AOMD::SPoint3(0,0,0),nullptr,0); // invalid parameter
}

int NullEdge::containsPoint(const AOMD::SPoint3 &pt) const
{  
  notImplemented();
  return 0;
}

int NullEdge::containsParam(double pt) const
{ 
  notImplemented();
  return 0;
}

GeoRep * NullEdge::geometry()
{ return new NullGeoRep(this,1); }

SVector3 NullEdge::firstDer(double par) const
{
  double der[3];
  notImplemented();
  return SVector3(der);
}

void NullEdge::nthDerivative(double param, int n, double *array) const
{
  notImplemented();
}

SPoint2 NullEdge::reparamOnFace(GFace * face, double epar, int dir) const
{
  notImplemented();
  return SPoint2(0,0);
}

double NullEdge::parFromPoint(const AOMD::SPoint3 &pt) const
{
  notImplemented();
  return MAXDOUBLE;
}


Logical::Value NullEdge::continuous(int) const
{ 
  return Logical::True;
}

Logical::Value NullEdge::degenerate(int) const
{ 
  return Logical::Unknown;
}

GeomType::Value NullEdge::geomType() const
{
  return GeomType::Unknown;
}

int NullEdge::geomDirection() const
{
  notImplemented();
  return 0;
}

double NullEdge::tolerance() const
{ return NullEntity::tolerance(); }


SSList<GEPoint> NullEdge::intersect(int fAxis, double fPar, GFace *f)
{
  SSList<GEPoint> intPts;
  notImplemented();  
  return intPts;
  
}

GVertex * NullEdge::split(double)
{
  notImplemented();
  return nullptr;
}

