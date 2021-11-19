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
#include "TPartialFace.h"
#include "SVector3.h"
#include <iostream>
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
#include "GRegion.h"
#include "GFPoint.h"
#include "GEPoint.h"
#include "TEdge.h"
#include "TFaceSplitEdge.h"
#include "SGModel.h"
#include "SSList.h"

using std::cerr;

TPartialFace::TPartialFace(SGModel *model, int t, GFace *f, int axis, 
			   int side, double par, GLoopUse *frontLoop, 
			   GShell *frontShell, GLoopUse *backLoop, 
			   GShell *backShell)
  : GFace(model,t,frontLoop,frontShell,backLoop,backShell), Axis(axis), 
    Side(side), Par(par), RealFace(f)
{}

TPartialFace::~TPartialFace()
= default;

AOMD::Range<double> TPartialFace::parBounds(int i) const
{ 
  AOMD::Range<double> b = RealFace->parBounds(i);
  if(i==Axis){ // if split on axis of request
    if(Side==0)
      b.high(Par);
    else
      b.low(Par);
  }
  return b;
}

int TPartialFace::paramDegeneracies(int dir, double *par)
  // degen. may not in in bounds of this part of face
{
  return RealFace->paramDegeneracies(dir,par);
}

AOMD::SBoundingBox3d TPartialFace::bounds() const
{ return RealFace->bounds(); }

SVector3 TPartialFace::normal(const SPoint2 &param) const
{ return RealFace->normal(param);}

Pair<SVector3,SVector3> TPartialFace::firstDer(const SPoint2 &param) const
{ return RealFace->firstDer(param); }

double * TPartialFace::nthDerivative(const SPoint2 &param, int n, double *array) const
{ return RealFace->nthDerivative(param,n,array); }

GFPoint TPartialFace::point(const SPoint2 &pt) const
{ return RealFace->point(pt); }

GFPoint TPartialFace::point(double par1,double par2) const
{ return RealFace->point(par1,par2); }

GFPoint TPartialFace::closestPoint(const AOMD::SPoint3 & queryPoint)
{ return RealFace->closestPoint(queryPoint); }

int TPartialFace::containsPoint(const AOMD::SPoint3 &pt) const
  // Note: this just returns what the real face returns
{
  return RealFace->containsPoint(pt);
}

int TPartialFace::containsParam(const SPoint2 &par) const
{
  if(RealFace->containsParam(par)){
    double p = par[Axis];
    if( (Side && p >= Par) || (!Side && p <= Par) )
      return 1;
  }
  return 0;
}

SPoint2 TPartialFace::parFromPoint(const AOMD::SPoint3 &pt) const
{ return RealFace->parFromPoint(pt); }

GeoRep * TPartialFace::geometry()
{ return RealFace->geometry(); }

Logical::Value TPartialFace::continuous(int dim) const
{ return Logical::True; }

Logical::Value TPartialFace::periodic(int dim) const
{ return Logical::False; }

Logical::Value TPartialFace::degenerate(int dim) const
{ return RealFace->degenerate(dim); }

GeomType::Value TPartialFace::geomType() const
{ return RealFace->geomType(); }

int TPartialFace::geomDirection() const
{ return RealFace->geomDirection(); }

double TPartialFace::tolerance() const
{ return RealFace->tolerance(); }

SSList<GEdge*> TPartialFace::splitU(double u)
{
  cerr << "TPartialFace::splitU - not done\n";
  return SSList<GEdge*>();
}

SSList<GEdge*> TPartialFace::splitV(double v)
{
  cerr << "TPartialFace::splitV - not done\n";
  return SSList<GEdge*>();
}

AOMD::GEntity * TPartialFace::getBaseRep()
{ return RealFace->getBaseRep(); }

Logical::Value TPartialFace::surfPeriodic(int dim) const
{ return RealFace->surfPeriodic(dim); }

//  Attribute * TPartialFace::attribute(const SString &type) const
//  {
//    if(RealFace)
//      return RealFace->attribute(type);
//    return 0;
//  }

//  SSList<Attribute *> TPartialFace::attributes(const SString &type) const
//  {
//    if(RealFace)
//      return RealFace->attributes(type);
//    return SSList<Attribute*>();
//  }

//  SSList<Attribute *> TPartialFace::attributes() const
//  {
//    if(RealFace)
//      return RealFace->attributes();
//    return SSList<Attribute*>();
//  }
