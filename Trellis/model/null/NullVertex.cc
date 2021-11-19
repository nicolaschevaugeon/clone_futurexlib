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
#include "NullVertex.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
//#include "GShell.h"
#include "GRegion.h"
#include "SSList.h"
#include "NullGeoRep.h"
#include "GVPoint.h"
#include "NullModel.h"

NullVertex::NullVertex(NullModel *m, int t) 
: GVertex(m,t), NullEntity() 
{}

NullVertex::~NullVertex() 
= default;

RepType::Value NullVertex::repType() const
{
  return RepType::Unknown;
}


GVPoint NullVertex::point() const
{
  double pt[3];
  //double param;
  notImplemented();
  return GVPoint(pt,this);
}

int NullVertex::containsPoint(const AOMD::SPoint3 &pt) const
{ 
  notImplemented();
  return 0;
}

GeoRep * NullVertex::geometry()
{ return new NullGeoRep(this,0); }

double NullVertex::tolerance() const
{ return NullEntity::tolerance(); }

