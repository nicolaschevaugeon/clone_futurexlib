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
#include "TEdgeSplitVertex.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
//#include "GShell.h"
#include "GRegion.h"
#include "SSList.h"
#include "GVPoint.h"
#include "GEPoint.h"
#include "GeoRep.h"

TEdgeSplitVertex::TEdgeSplitVertex(SGModel *model, int t, GEdge *e, double param)
: GVertex(model,t), SplitEdge(e), Param(param)
{}

TEdgeSplitVertex::~TEdgeSplitVertex() 
= default;

GVPoint TEdgeSplitVertex::point() const
{ return GVPoint(SplitEdge->point(Param),this); }

int TEdgeSplitVertex::containsPoint(const AOMD::SPoint3 &pt) const
{ 
  AOMD::SPoint3 vpt = SplitEdge->point(Param);
  return (dist(vpt,pt) < tolerance());
}

GeoRep * TEdgeSplitVertex::geometry()
{ 
  return nullptr;
  //return new GeoRep(this,0); 
}

double TEdgeSplitVertex::tolerance() const
{ return SplitEdge->tolerance(); }

//  Attribute * TEdgeSplitVertex::attribute(const SString &type) const
//  {
//    if(SplitEdge)
//      return SplitEdge->attribute(type);
//    return 0;
//  }

//  SSList<Attribute *> TEdgeSplitVertex::attributes(const SString &type) const
//  {
//    if(SplitEdge)
//      return SplitEdge->attributes(type);
//    return SSList<Attribute*>();
//  }

//  SSList<Attribute *> TEdgeSplitVertex::attributes() const
//  {
//    if(SplitEdge)
//      return SplitEdge->attributes();
//    return SSList<Attribute*>();
//  }
