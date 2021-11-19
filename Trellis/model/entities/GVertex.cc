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
#include "GVertex.h"
#include "GEdgeUsePair.h"
#include "GVertexUse.h"
#include "GLoopUse.h"
#include "GFaceUse.h"
#include "GFace.h"
#include "GVPoint.h"
#include "SGModel.h"
#include "MessageOut.h"
#include "EDList.h"
#include "GEdge.h"
#include "SString.h"

using std::ostream;
using std::istream;
using std::cerr;

GVertex::GVertex(SGModel *m, int t) 
: GEntity(m,t)
{ }

TopoType::Value GVertex::type() const
{ return TopoType::Vertex; }

int GVertex::dim() const
{ return 0; }

AOMD::SBoundingBox3d GVertex::bounds() const
{ return AOMD::SBoundingBox3d(point()); }

SSList<GRegion*> GVertex::regions() const
{
  SSList<GRegion*> rs;
  GVertexUse *vu;
  for(SSListCIter<GVertexUse*> vuIter(Uses); vuIter(vu);){
    GEdgeUsePair *eu;
    for(SSListCIter<GEdgeUsePair*> euIter = vu->firstEdgeUse(); euIter(eu); ){
      rs.appendUnique(eu->loopUse(0)->faceUse()->face()->regions());
      rs.appendUnique(eu->loopUse(1)->faceUse()->face()->regions());
    }
  }
  return rs;
}

SSList<GFace*> GVertex::faces() const
{
  SSList<GFace*> fs;
  GVertexUse *vu;
  for(SSListCIter<GVertexUse*> vuIter(Uses); vuIter(vu);){
    GEdgeUsePair *eu;
    for(SSListCIter<GEdgeUsePair*> euIter = vu->firstEdgeUse(); euIter(eu); ){
      fs.appendUnique(eu->loopUse(0)->faceUse()->face());
      fs.appendUnique(eu->loopUse(1)->faceUse()->face());
    }
  }
  return fs;
}

SSList<GEdge*> GVertex::edges() const
{
  SSList<GEdge*> es;
  GVertexUse *vu;
  for(SSListCIter<GVertexUse*> vuIter(Uses); vuIter(vu);)
    es.appendUnique(vu->edges());
  return es; 
}

SSList<GVertex*> GVertex::vertices() const
{
  SSList<GVertex*> vs;
  return vs;
}

int GVertex::classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const
{
  if(containsPoint(pt))
    return 1;
  return 0;
}

int GVertex::containsPoint(const AOMD::SPoint3 &pt) const
{
  AOMD::SPoint3 vpt = point();
  double distsq = (pt[0]-vpt[0])*(pt[0]-vpt[0])+(pt[1]-vpt[1])*(pt[1]-vpt[1])+
              (pt[2]-vpt[2])*(pt[2]-vpt[2]);
  double tol = tolerance();
  if(distsq <= tol*tol)
    return 1;
  return 0;
}

GVertexUse * GVertex::addEdge(GEdgeUsePair *e)
{ 
  GVertexUse *vu = new GVertexUse(model(),this,e);
  Uses.append(vu); 
  return vu;
}

int GVertex::numEdges() const
{ 
  return edges().size(); 
}

int GVertex::inClosure(GEntity *ent) const
{
  return this==ent;
}

void GVertex::removeUse(GVertexUse *vu)
{
  Uses.remove(vu);
}

int GVertex::numUses() const
{ return Uses.size(); }

SSList<GVertexUse*> GVertex::uses() const
{ return Uses; }

void GVertex::write(ostream &out)
{
  int i = 0;
  GVertexUse *vu;
  for(SSListCIter<GVertexUse*> vuIter(Uses); vuIter(vu); )
    vu->setID(i++);
  out << tag() << " " << id() << " " << Uses.size() << "\n";
}

void GVertex::read(istream &in)
{
  int nuse, id;
  in >> id >> nuse;
  setID(id);
  for(int i = 0; i < nuse; i++){
    GVertexUse *vu = new GVertexUse(model(),this);
    Uses.append(vu);
    vu->setID(i);
  }
}

GVertexUse * GVertex::use(int id)
{ return Uses.nth(id); }

void GVertex::info() const
{
  cerr << name() << " - " << (void*)this << "\n";
  GVertex *vtx = (GVertex*)this;
  GVPoint point = vtx->point();
  cerr << "tag = " << vtx->tag() << "\n";
  cerr << "  number of uses = " << vtx->numUses() << "\n";
  cerr << "  " << point.x() << ", " << point.y() << ", " << point.z() << "\n";
  
  SSList<GEdge*> ve = vtx->edges();
  cerr << "  " << ve.size() << " edges:  ";
  GEdge *e;
  for(SSListCIter<GEdge*> eIter(ve); eIter(e); ){
    cerr << e->tag() << " ";
  }
  cerr << "\n";
  SSList<GFace*> vf = vtx->faces();
  cerr << "  " << vf.size() << " faces\n";
  SSList<GRegion*> vr = vtx->regions();
  cerr << "  " << vr.size() << " regions\n";
}






