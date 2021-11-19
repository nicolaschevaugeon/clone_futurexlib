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
#include "TVertex.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
//#include "GShell.h"
#include "GRegion.h"
#include "SSList.h"
#include "GVPoint.h"

TVertex::TVertex(SGModel *model, int t, GVertex *v)
: GVertex(model,t), RealVertex(v)
{}

TVertex::~TVertex() 
= default;

GVPoint TVertex::point() const
{ return RealVertex->point(); }

GeoRep * TVertex::geometry()
{ return RealVertex->geometry(); }

double TVertex::tolerance() const
{ 
  if (RealVertex)
    return RealVertex->tolerance(); 
  return 1.0e-6;
}

int TVertex::containsPoint(const AOMD::SPoint3 &pt) const
{ return RealVertex->containsPoint(pt); }

//  Attribute * TVertex::attribute(const SString &type) const
//  {
//    if(RealVertex)
//      return RealVertex->attribute(type);
//    return  SModelMember::attribute(type); 
//    return 0;
//  }

//  SSList<Attribute *> TVertex::attributes(const SString &type) const
//  {
//    if(RealVertex)
//      return RealVertex->attributes(type);
//    return  SModelMember::attributes(type); 
//  }

//  SSList<Attribute *> TVertex::attributes() const
//  {
//    if(RealVertex)
//      return RealVertex->attributes();
//    return  SModelMember::attributes(); 
//  }
