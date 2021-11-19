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
#include "GLoopUse.h"
#include "GEdge.h"
#include "GEdgeUsePair.h"
#include "GEdgeUse.h"
#include "GVertexUse.h"
#include "SGModel.h"
#include "SSList.h"
using std::cerr;
using std::ostream;
using std::istream;


GLoopUse::GLoopUse(SGModel *m, GFaceUse *fu, const SSList<GEdge*> &edges, const SSList<int> &dirs)
: GEntity(m,0), FaceUse(fu)
{
  int dir;
  GEdge *edge;
  SSListCIter<GEdge *> eIter(edges);
  SSListCIter<int> dIter(dirs);
  for( ; eIter(edge), dIter(dir) ; )
    appendEdge(edge,dir);
  mergeVertexUses();
}

GLoopUse::GLoopUse(SGModel *m, const SSList<GEdgeUse*> & edges, 
		   GFaceUse *faceUse, GLoopUse *old)
  : GEntity(m,0), Edges(edges), FaceUse(faceUse)
{
  GEdgeUse *eu;
  for( SSListIter<GEdgeUse*> euIter(Edges); euIter(eu); )
    eu->setLoopUse(this);
}

GLoopUse::GLoopUse(SGModel *m, GFaceUse *faceUse)
  : GEntity(m,0)
{
}

GLoopUse::~GLoopUse()
= default;

void GLoopUse::setMate(GLoopUse *m)
  // set the loop use mate for this use, this must be called by the face
  // that creates the loop uses
{
  d_mate = m;
}

GeoRep * GLoopUse::geometry()
{ 
  cerr << "GLoopUse::geometry() called\n"; 
  return nullptr;
}

int GLoopUse::numEdges() const
{ return Edges.size(); }


SSListCIter<GEdgeUse*> GLoopUse::firstEdgeUse() const
{
  return SSListCIter<GEdgeUse*>(Edges);
}

SSList<GEdge*> GLoopUse::edges() const
{ 
  SSList<GEdge*> eds;
  GEdgeUse *eu;
  for(SSListCIter<GEdgeUse*> eIter(Edges); eIter(eu); )
    eds.append(eu->edge());
  return eds; 
}

GFaceUse* GLoopUse::faceUse() const
{ return FaceUse; }

//void GLoopUse::setMesh(Mesh *, MEntity *ent)
//{ cerr << "GLoopUse::setMesh() called\n"; }

TopoType::Value GLoopUse::type() const
{ return TopoType::LoopUse; }

int GLoopUse::dim() const
{ return 2; }

AOMD::SBoundingBox3d GLoopUse::bounds() const
// for a loop bounding box is bounds of all edges
{
  AOMD::SBoundingBox3d box;
  GEdge *e;
  for(SSListCIter<GEdge*> eIter(edges()); eIter(e); )
    box += e->bounds();
  return box;
}

void GLoopUse::appendEdge(GEdge *e, int dir)
{
  //int which = Edges.size();
  Edges.append(e->addLoop(this,dir));  // returns a GEdgeUse for e
}

void GLoopUse::setFaceUse(GFaceUse *f)
{
  FaceUse = f;
}

void GLoopUse::mergeVertexUses()
{
  // make sure case of only one edge in loop is handled correctly
  GEdgeUse *eus, *eusLast;
  //int dir,lastDir;

  eusLast = Edges.nth(Edges.size()-1);

  for(SSListCIter<GEdgeUse*> euIter(Edges); euIter(eus); ){
    GVertexUse *vu0 = eusLast->vertexUse();
    GVertexUse *vu1 = eus->otherSide()->vertexUse();
    if(vu1 != vu0)
      vu1->mergeWith(vu0);
    eusLast = eus;
  }
}

int GLoopUse::inClosure(GEntity *ent) const
{
  if(this==ent) 
    return 1;
  GEdgeUse *e;
  for(SSListCIter<GEdgeUse*> eIter(Edges); eIter(e); )
    if(e->inClosure(ent))
      return 1;
  return 0;
}

void GLoopUse::addBefore(GEdgeUse *eu, GEdgeUse *newUse)
{
  GEdgeUse *e;
  for(SSListIter<GEdgeUse*> euIter(Edges); euIter(e); ){
    if(e==eu)
      euIter.insertBefore(newUse);
  }
}

void GLoopUse::addAfter(GEdgeUse *eu, GEdgeUse *newUse)
{
  GEdgeUse *e;
  for(SSListIter<GEdgeUse*> euIter(Edges); euIter(e); ){
    if(e==eu)
      euIter.insert(newUse);
  }
}

void GLoopUse::replace(GEdgeUse *eu, GEdgeUse *newUse)
{
  GEdgeUse *e;
  for(SSListIter<GEdgeUse*> euIter(Edges); euIter(e); ){
    if(e==eu)
      euIter.replaceCurrent(newUse);
  }
}

GLoopUse * GLoopUse::split(GEdge *e, int dir, GEdgeUse *after0, GEdgeUse *after1)
  // split this loop use with the introduction of new uses of edge e
  // into each half of the loopuse. 
  // the first vertex of e is after after0 and the second is after after1
  // after0 is before after1 in the list of edgeUseSides
{
  SSList<GEdgeUse*> otherLoopEdges;

  GEdgeUse *eus;
  SSListIter<GEdgeUse*> eusIter(Edges);
  while(eusIter(eus))
    if(eus==after0)
      break;
  while(eusIter(eus)){
    otherLoopEdges.append(eus);
    eusIter.remove();
    if(eus==after1)
      break;
  }
  GEdgeUse *newUse = e->addLoop(this,dir);
  Edges.append(newUse); 

  GLoopUse *newLoopUse = new GLoopUse(model(),otherLoopEdges,FaceUse,this);
  newLoopUse->Edges.append(newUse->otherSide());
  newUse->otherSide()->setLoopUse(newLoopUse);

  return newLoopUse;
}

void GLoopUse::merge(GEdge *e, int dir, GEdgeUse *after0, GEdgeUse *after1)
  // merge two loops by the introduction of edge e. after0 is the use in this
  // after which the edge is inserted. after1 is the use in the other loop
  // after which the edge is inserted.
{
  GLoopUse *otherLoop = after1->loopUse();
#ifdef DEBUG
  if(this == otherLoop)
    cerr << "GLoopUse::merge - same loop\n";
  if(faceUse() != otherLoop->faceUse())
    cerr << "GLoopUse::merge - different face uses\n";
#endif
  // first, find the location in each list of uses for the insertion
  GEdgeUse *eus1,*eus2;
  int found1=0, found2=0;
  SSListIter<GEdgeUse*> eu1Iter(Edges);
  while(eu1Iter(eus1))
    if(eus1 == after0){
      found1 = 1;
      break;
    }
  SSListIter<GEdgeUse*> eu2Iter(otherLoop->Edges);
  while(eu2Iter(eus2))
    if(eus2 == after1){
      found2 = 1;
      break;
    }
  if(! (found1 &&found2))
    cerr << "GLoopUse::merge - error 1\n";
  // NOTE: need to take care of new vertex uses

  // add uses to this loop and change ownership
  GEdgeUse *newEus = e->addLoop(this,dir);
  eu1Iter.insert(newEus);
  eu1Iter(); // advance iterator
  while(eu2Iter(eus2)){
    eu1Iter.insert(eus2);
    eu1Iter(); // make sure to advance the iterator
    eus2->setLoopUse(this);
  }
  eu2Iter.reset();
  do{
    eu2Iter(eus2);
    eu1Iter.insert(eus2);
    eu1Iter(); // advance the iterator
    eus2->setLoopUse(this);
  } while(eus2!=after1);
  eu1Iter.insert(newEus->otherSide());
  newEus->otherSide()->setLoopUse(this);
  // now this loopuse should have everything in it
}

void GLoopUse::write(ostream &out)
{
  GEdgeUse *eus;
  out << Edges.size() << " ";
  for(SSListCIter<GEdgeUse*> eusIter(Edges); eusIter(eus); ){
    out << eus->edge()->id() << " " << eus->use()->id() << " ";
    out << eus->dir() << " ";
  }
}

void GLoopUse::read(istream &in)
{
  int nedge;
  in >> nedge;
  for(int i = 0; i < nedge; i++){
    int edgeid, useid, dir;
    in >> edgeid >> useid >> dir;
    GEdgeUse *eus = model()->edgeByID(edgeid)->use(useid)->side(dir);
    Edges.append(eus);
    eus->setLoopUse(this);
  }
}


void GLoopUse::info() const
{
  cerr << "loopuse: ";
  SSListCIter<GEdgeUse*> eusIter(Edges);
  GEdgeUse * eus;
  while( eusIter(eus) ){
    cerr << "(" << eus->dir() << "," << eus->edge()->tag() << ") ";
  }
  cerr << "\n";

}
