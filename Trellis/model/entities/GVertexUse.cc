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
#include "GVertexUse.h"
#include "GLoopUse.h"
#include "GVertex.h"
#include "GEdgeUsePair.h"
using std::cerr;
GVertexUse::GVertexUse(SGModel *m, GVertex *v, GEdgeUsePair *eu)
  : GEntity(m,0), OwningVertex(v)
{
  EdgeUses.append(eu);
}

GVertexUse::GVertexUse(SGModel *m, GVertex *v)
  : GEntity(m,0), OwningVertex(v)
{
}

GVertexUse::~GVertexUse()
= default;

GeoRep * GVertexUse::geometry()
{ 
  cerr << "GVertexUse::geometry() called\n"; 
  return nullptr;
}

AOMD::SBoundingBox3d GVertexUse::bounds() const
{
  return OwningVertex->bounds();
}

int GVertexUse::inClosure(GEntity *ent) const
{
  if(this==ent || OwningVertex==ent)
    return 1;
  return 0;
}

TopoType::Value GVertexUse::type() const
{ return TopoType::VertexUse; }

int GVertexUse::dim() const
{ return 0; }

SSList<GEdge*> GVertexUse::edges() const
{
  SSList<GEdge*> es;
  GEdgeUsePair *eu;
  for(SSListCIter<GEdgeUsePair*> euIter(EdgeUses); euIter(eu);)
    es.appendUnique(eu->edge());
  return es; 
}

GShell * GVertexUse::shell() const
{
  GShell *su =nullptr;
  if(EdgeUses.size()){
    su = EdgeUses.nth(0)->shell();
  }
  return su;
}

SSList<GEdgeUsePair*> GVertexUse::edgeUses() const
{ return EdgeUses; }

SSList<GFaceUse*> GVertexUse::faceUses() const
{
  SSList<GFaceUse*> fu;
  GEdgeUsePair *eu;
  for(SSListCIter<GEdgeUsePair*> euIter(EdgeUses); euIter(eu); ){
    fu.appendUnique(eu->loopUse(0)->faceUse());
    fu.appendUnique(eu->loopUse(1)->faceUse());
  }
  return fu;
}


GVertex * GVertexUse::vertex() const
{ return OwningVertex; }

SSListCIter<GEdgeUsePair*> GVertexUse::firstEdgeUse() const
{ return SSListCIter<GEdgeUsePair*>(EdgeUses); }

void GVertexUse::mergeWith(GVertexUse *vu)
{ 
  if(this!=vu){
    if(OwningVertex != vu->vertex())
      cerr << "GVertexUse::mergeWith - error\n";
    OwningVertex->removeUse(vu);
    GEdgeUsePair *eu;
    for(SSListIter<GEdgeUsePair*> euIter(vu->EdgeUses); euIter(eu); ){
      EdgeUses.append(eu);
      if(eu->vertexUse(0)==vu)
	eu->setVertexUse(0,this);
      if(eu->vertexUse(1)==vu)
	eu->setVertexUse(1,this);
    }
  }
}

void GVertexUse::replaceEdgeUse(GEdgeUsePair *old, GEdgeUsePair *newUse)
  // replace the first occurance of the use "old" with the use "new"
{
  GEdgeUsePair *eu;
  for(SSListIter<GEdgeUsePair*> euIter(EdgeUses); euIter(eu); ){
    if(eu==old)
      euIter.replaceCurrent(newUse);
  }
}

void GVertexUse::addEdgeUse(GEdgeUsePair *eu)
{ EdgeUses.append(eu); }

void GVertexUse::removeEdgeUse(GEdgeUsePair *eu)
{ EdgeUses.remove(eu); }
