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
#include "GEntity.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
//#include "GShell.h"
#include "GRegion.h"
//#include "SSList.h"
#include "SGModel.h"
#include "MessageOut.h"
#include <cstdio>

namespace AOMD {

GEntity::GEntity(SGModel *model, int t)
  : SModelMember(model,t)
{}

GEntity::~GEntity() 
= default;

SString GEntity::name() const
{
  SString n;
  switch(dim()){
  case 0:
    n = "vertex ";
    break;
  case 1:
    n = "edge ";
    break;
  case 2:
    n = "face ";
    break;
  case 3:
    n = "region ";
    break;
  default:
    n = "unknown ";
    break;
  }
  char tagstr[20];
  sprintf(tagstr,"%d",tag());
  return n+SString(tagstr);
}

RepType::Value GEntity::repType() const
{
  return RepType::Parametric;
}

double GEntity::tolerance() const
{ 
  notImplemented();
  return -1; 
}

Logical::Value GEntity::continuous(int ) const
{ 
  notImplemented();
  return Logical::Unknown; 
}

Logical::Value GEntity::periodic(int ) const
{ 
  notImplemented();
  return Logical::Unknown; 
}

Logical::Value GEntity::degenerate(int ) const
{ 
  notImplemented();
  return Logical::Unknown; 
}

GeomType::Value GEntity::geomType() const
{ return GeomType::Unknown; }

Range<double> GEntity::parBounds(int ) const
{ return Range<double> (0,0); }

int GEntity::containsPoint(const SPoint3 &pt) const
{ 
  InternalError("GEntity::containsPoint called");
  return 0;
}

int GEntity::classifyPoint(const SPoint3 &pt, int chkb, GEntity **bent) const
{ 
  InternalError("GEntity::classifyPoint called");
  return 0;
}

int GEntity::geomDirection() const
{ return 1; }

SSList<GRegion*> GEntity::regions() const
{
  notImplemented();
  return SSList<GRegion*>();
}

SSList<GFace*> GEntity::faces() const
{
  notImplemented();
  return SSList<GFace*>();
}

SSList<GEdge*> GEntity::edges() const
{
  notImplemented();
  return SSList<GEdge*>();
}

SSList<GVertex*> GEntity::vertices() const
{
  notImplemented();
  return SSList<GVertex*>();
}

GEntity * GEntity::getBaseRep()
{ return this; }

void * GEntity::getNativePtr() const
{ 
  return nullptr;
}

int GEntity::getNativeInt() const
{ 
  InternalError("GEntity::getNativeInt called");
  return 0;
}

}
