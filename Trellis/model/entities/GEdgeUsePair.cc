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
#include "GEdgeUsePair.h"
#include "GEdge.h"
#include "GVertex.h"
#include "GVertexUse.h"
#include "GFaceUse.h"
#include "GLoopUse.h"
#include "SGModel.h"

using std::cerr;
using std::ostream;
using std::istream;



GEdgeUsePair::GEdgeUsePair(SGModel *m, GEdge *e, GLoopUse *l, int dir, GVertex *v0, GVertex *v1)
  : GEntity(m,0), OwningEdge(e), PosSide(m,this), NegSide(m,this)
  // edge uses are generally constructed with only one loop attached and
  // one side used
{
  dir ? PosSide.setLoopUse(l) : NegSide.setLoopUse(l);
  if(v0)
    NegSide.setVertexUse( v0->addEdge(this) );
  if(v1)
    PosSide.setVertexUse( v1->addEdge(this) );
}

GEdgeUsePair::GEdgeUsePair(SGModel *m, GEdge *e)
  : GEntity(m,0), OwningEdge(e), PosSide(m,this,nullptr,nullptr), NegSide(m,this,nullptr,nullptr)
{ }

GEdgeUsePair::~GEdgeUsePair()
= default;

GeoRep * GEdgeUsePair::geometry()
{ 
  cerr << "GEdgeUsePair::geometry() called\n"; 
  return nullptr;
}

AOMD::SBoundingBox3d GEdgeUsePair::bounds() const
{
  return OwningEdge->bounds();
}

int GEdgeUsePair::inClosure(GEntity *ent) const
  // not sure what this should do
{
  if(this==ent || OwningEdge==ent || vertexUse(0)->inClosure(ent) ||
     vertexUse(1)->inClosure(ent))
    return 1;
  return 0;
}

TopoType::Value GEdgeUsePair::type() const
{ return TopoType::EdgeUse; }

int GEdgeUsePair::dim() const
{ return 1; }

GEdge * GEdgeUsePair::edge() const
{ return OwningEdge; }

GShell * GEdgeUsePair::shell() const
{
  GShell *su =nullptr;
  su = NegSide.shell();
  if(!su)
    su = PosSide.shell();
  return su;
}

GLoopUse * GEdgeUsePair::loopUse(int dir) const
{ return dir ? PosSide.loopUse() : NegSide.loopUse(); }

GVertexUse * GEdgeUsePair::vertexUse(int dir) const
{ return dir ? PosSide.vertexUse() : NegSide.vertexUse(); }

GEdgeUse * GEdgeUsePair::side(int which)
{ return which ? &PosSide : &NegSide; }

int GEdgeUsePair::whichSide(const GEdgeUse *eus) const
{ 
  if(eus == &PosSide)
    return 1;
  if(eus == &NegSide)
    return 0;
  cerr << "GEdgeUsePair::whichSide - error\n";
  return -1;
}

GEdgeUse * GEdgeUsePair::otherSide(GEdgeUse *eus)
{ 
  if(eus == &PosSide)
    return &NegSide;
  if(eus == &NegSide)
    return &PosSide;
  cerr << "GEdgeUsePair::otherSide - error\n";
  return nullptr;
}


void GEdgeUsePair::setLoopUse(GLoopUse *lu, int dir)
{ dir ? PosSide.setLoopUse(lu) : NegSide.setLoopUse(lu); }

void GEdgeUsePair::setVertexUse(int dir, GVertexUse *vu)
{
  if(dir)
    PosSide.setVertexUse(vu);
  else
    NegSide.setVertexUse(vu);
}

void GEdgeUsePair::mergeWith(GEdgeUsePair *eu)
  // merge two uses, update loopuses
{ 
  GLoopUse *lu{nullptr};
  GEdgeUse *oldeu, *neweu;

  if(OwningEdge != eu->edge())
    cerr << "GEdgeUsePair::mergeWith - error\n";
  OwningEdge->removeUse(eu);
  if(loopUse(0) && !loopUse(1) && eu->loopUse(1) && !(eu->loopUse(0))){
    neweu = &PosSide;
    oldeu = eu->side(1);
    lu = eu->loopUse(1);
    PosSide.setLoopUse(lu);
  } else if(loopUse(1) && !loopUse(0) && eu->loopUse(0) && !(eu->loopUse(1))){
    neweu = &NegSide;
    oldeu = eu->side(0);
    lu = eu->loopUse(0);
    NegSide.setLoopUse(lu);
  }  else {
    cerr << "GEdgeUsePair::mergeWith - problem with uses\n";
    cerr << "Problem edge tag : " << edge()->tag() << "\n";
    cerr << "Tags of faces of uses being merged : ";
    if(loopUse(0))
      cerr << loopUse(0)->faceUse()->face()->tag() << " ";
    if(loopUse(1))
      cerr << loopUse(1)->faceUse()->face()->tag() << " ";
    if(eu->loopUse(0))
      cerr << eu->loopUse(0)->faceUse()->face()->tag() << " ";
    if(eu->loopUse(1))
      cerr << eu->loopUse(1)->faceUse()->face()->tag() << " ";
    cerr << "\n";
  }
  eu->vertexUse(0)->removeEdgeUse(eu);
  eu->vertexUse(1)->removeEdgeUse(eu);
  vertexUse(0)->mergeWith(eu->vertexUse(0));
  vertexUse(1)->mergeWith(eu->vertexUse(1));
  lu->replace(oldeu,neweu);
}

GEdgeUsePair * GEdgeUsePair::split(GVertex *v, GEdge *se)
  // split edge use since the underlying edge was just split
  // the new vertex and new edge are passed in, the new half
  // of the edge use is on the second half of the edge
{
  GEdgeUsePair * eu = new GEdgeUsePair(model(),se);
  //GVertexUse *vu0 = NegSide.vertexUse();
  GVertexUse *vu1 = PosSide.vertexUse();
  vu1->replaceEdgeUse(this,eu);
  eu->PosSide.setVertexUse(vu1);
  GVertexUse *vuc = v->addEdge(this); // creates new vertex use
  PosSide.setVertexUse(vuc);
  vuc->addEdgeUse(eu);
  eu->NegSide.setVertexUse(vuc);

  // now need to update the loop uses
  if(loopUse(0) == loopUse(1))
    cerr << "Warning: GEdgeUsePair::split - check this case\n";
  loopUse(0)->addBefore(&NegSide,&(eu->NegSide));
  eu->setLoopUse(loopUse(0),0);
  loopUse(1)->addAfter(&PosSide,&(eu->PosSide));
  eu->setLoopUse(loopUse(1),1);

  return eu;
}

void GEdgeUsePair::replaceOwner(GEdge *enew)
{
  OwningEdge = enew;
}

void GEdgeUsePair::write(ostream &out)
{
  for(int i=0; i < 2; i++){
    GVertexUse *vu = vertexUse(i);
    out << vu->vertex()->id() << " " << vu->id() << " ";
  }
}

void GEdgeUsePair::read(istream &in)
{
  int vid, vuid;
  in >> vid >> vuid;
  GVertexUse *vu = model()->vertexByID(vid)->use(vuid);
  NegSide.setVertexUse(vu);
  vu->addEdgeUse(this);
  in >> vid >> vuid;
  vu = model()->vertexByID(vid)->use(vuid);
  PosSide.setVertexUse(vu);
  vu->addEdgeUse(this);
}

void GEdgeUsePair::info() const
{
  cerr << "  " << (void*)this << "\n";
  PosSide.info();
  NegSide.info();
}
