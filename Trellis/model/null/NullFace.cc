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
#include "NullFace.h"
#include "NullEdge.h"
#include "NullRegion.h"
#include "SVector3.h"
#include <iostream>
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
//#include "GShell.h"
#include "GRegion.h"
#include "GEPoint.h"
#include "GFPoint.h"
#include "SSList.h"
#include "NullGeoRep.h"
#include "NullModel.h"

#include <cfloat>
#define MAXDOUBLE DBL_MAX

NullFace::NullFace(NullModel *m, int t, const SSList<GEdge*> &edges, const SSList<int> &dirs) 
  : GFace(m,t, edges,dirs), NullEntity() 
{}

NullFace::NullFace(NullModel *m, int t) 
  : GFace(m,t), NullEntity() 
{}

NullFace::~NullFace()
= default;

RepType::Value NullFace::repType() const
{
  return RepType::Unknown;
}

AOMD::Range<double> NullFace::parBounds(int i) const
{ 
  notImplemented();
  return( AOMD::Range<double>(0,0) );
}


int NullFace::paramDegeneracies(int dim, double *par)
{

  notImplemented();
  return 0;
}

AOMD::SBoundingBox3d NullFace::bounds() const
{
  AOMD::SBoundingBox3d box;
  notImplemented();
  return box;
}

SVector3 NullFace::normal(const SPoint2 &param) const
{
  double norm[3];
  notImplemented();
  return SVector3(norm);
}

Pair<SVector3,SVector3> NullFace::firstDer(const SPoint2 &param) const
{
  double tanSpace[6];
  notImplemented();
  return Pair<SVector3,SVector3>( SVector3(tanSpace),SVector3(&(tanSpace[3])));
}

double * NullFace::nthDerivative(const SPoint2 &param, int n, double *array) const
{
  notImplemented();
  return nullptr;
}

SPoint2 NullFace::parFromEdgePar(GEdge *edge, double epar, int dir) const
{
  notImplemented();
  return SPoint2(0,0);
}

GFPoint NullFace::point(const SPoint2 &pt) const
{ return point(pt.x(),pt.y()); }

GFPoint NullFace::point(double par1,double par2) const
{
  notImplemented();
  return GFPoint(AOMD::SPoint3(0,0,0),nullptr,SPoint2(par1,par2)); // invalid parameter
}

GFPoint NullFace::closestPoint(const AOMD::SPoint3 & queryPoint)
{
  notImplemented();
  return GFPoint(AOMD::SPoint3(0,0,0), this, SPoint2(0,0));
}

int NullFace::containsPoint(const AOMD::SPoint3 &pt) const
{ 
  notImplemented();
  return 0;
}

int NullFace::containsParam(const SPoint2 &pt) const
{ 
  notImplemented();
  return 0;
}

SPoint2 NullFace::parFromPoint(const AOMD::SPoint3 &pt) const
{
  notImplemented();
  return SPoint2(MAXDOUBLE,MAXDOUBLE);
}


GeoRep * NullFace::geometry()
{ return new NullGeoRep(this,2); }

Logical::Value NullFace::continuous(int) const
{ 
  return Logical::True;
}

Logical::Value NullFace::periodic(int dim) const
{ 
  return Logical::True;
}

Logical::Value NullFace::degenerate(int) const
{ 
  return Logical::Unknown;
}

GeomType::Value NullFace::geomType() const
{
  return GeomType::Unknown;
}

int NullFace::geomDirection() const
{
    notImplemented();
    return 1;
}

double NullFace::tolerance() const
{ return NullEntity::tolerance(); }

SVector3 NullFace::inwardSecant(GEdge *e, double ept, int dir)
{
  notImplemented();
  return SVector3(0,0,0);
}

Logical::Value NullFace::surfPeriodic(int dim) const
{ 
  notImplemented();
  return Logical::False; 
}
