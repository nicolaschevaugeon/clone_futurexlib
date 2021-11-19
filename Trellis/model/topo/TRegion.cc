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
#include "TRegion.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
#include "GShell.h"
#include "GRegion.h"

#include "SSList.h"

TRegion::TRegion(SGModel *model, int t, GRegion *r,
		 const SSList<GFace*> & faces, const SSList<int> &dirs)
: GRegion(model, t, faces,dirs), RealRegion(r)
{}

TRegion::TRegion(SGModel *model, int t)
: GRegion(model,t), RealRegion(nullptr)
{}

TRegion::~TRegion() 
= default;


GeoRep * TRegion::geometry()
{ return RealRegion->geometry(); }

double TRegion::tolerance() const
{ 
  if (RealRegion)
    return RealRegion->tolerance(); 
  return 1.0e-6;

}

int TRegion::containsPoint(const AOMD::SPoint3 &pt) const
{ return RealRegion->containsPoint(pt); }

//  Attribute * TRegion::attribute(const SString &type) const
//  {
//    if(RealRegion)
//      return RealRegion->attribute(type);
//    return  SModelMember::attribute(type); 
//  }

//  SSList<Attribute *> TRegion::attributes(const SString &type) const
//  {
//    if(RealRegion)
//      return RealRegion->attributes(type);
//    return  SModelMember::attributes(type); 
//  }

//  SSList<Attribute *> TRegion::attributes() const
//  {
//    if(RealRegion)
//      return RealRegion->attributes();
//    return  SModelMember::attributes(); 
//  }
