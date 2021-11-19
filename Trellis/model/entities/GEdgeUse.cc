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
#include "GEdgeUse.h"
#include "GEdge.h"
#include "GFace.h"
#include "GVertex.h"
#include "GVertexUse.h"
#include "GFaceUse.h"
#include "GLoopUse.h"
#include "GEdgeUsePair.h"

#include "SString.h"
using std::cerr;

GEdgeUse::GEdgeUse(SGModel *m, GEdgeUsePair *ep, GLoopUse *l, GVertex *v0)
  : GEntity(m,0), Use(ep)
{
  LoopUse = l;
  if(v0)
    VertexUse = v0->addEdge(ep);
}

GEdgeUse::GEdgeUse(SGModel *m, GEdgeUsePair *e)
  : GEntity(m,0), Use(e)
{
  LoopUse = nullptr;
  VertexUse = nullptr;
}

GEdgeUse::~GEdgeUse()
= default;

GeoRep * GEdgeUse::geometry()
{ 
  cerr << "GEdgeUse::geometry() called\n"; 
  return nullptr;
}

AOMD::SBoundingBox3d GEdgeUse::bounds() const
{
  return Use->bounds();
}

int GEdgeUse::inClosure(GEntity *ent) const
{ return Use->inClosure(ent); }

TopoType::Value GEdgeUse::type() const
{ 
  notImplemented();
  return TopoType::EdgeUse; 
}

int GEdgeUse::dim() const
{ return 1; }

GEdge * GEdgeUse::edge() const
{ return Use->edge(); }

int GEdgeUse::dir() const
{ return Use->whichSide(this); }

GEdgeUsePair * GEdgeUse::use() const
{ return Use; }

GEdgeUse * GEdgeUse::otherSide()
{ return Use->otherSide(this); }

GShell * GEdgeUse::shell() const
{
  GShell *su =nullptr;
  if(LoopUse){
    GFaceUse *fu = LoopUse->faceUse();
    if(fu)
      su = fu->shell();
  }
  return su;
}

GLoopUse * GEdgeUse::loopUse() const
{ return LoopUse; }

GVertexUse * GEdgeUse::vertexUse() const
{ return VertexUse; }

void GEdgeUse::setLoopUse(GLoopUse *lu)
{ LoopUse = lu; }

void GEdgeUse::setVertexUse(GVertexUse *vu)
{ VertexUse = vu; }


void GEdgeUse::info() const
{
  cerr << "    " << (void*)this << "\n";
  cerr << "    loop use: " << LoopUse << " ->face use: " << LoopUse->faceUse();
  cerr << " ->face: " << LoopUse->faceUse()->face();
  cerr << "(" << LoopUse->faceUse()->face()->name() << ") dir= ";
  cerr << LoopUse->faceUse()->dir() << "\n";
}
