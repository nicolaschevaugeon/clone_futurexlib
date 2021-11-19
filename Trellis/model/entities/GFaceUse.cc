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
#include "GFaceUse.h"
#include "GLoopUse.h"
#include "GFace.h"
#include "GEdge.h"
#include "GEdgeUse.h"
#include "GVertexUse.h"

using std::cerr;
using std::ostream;
using std::istream;

GFaceUse::GFaceUse(SGModel *m, GFace * face)
  : GEntity(m,0), Face(face), Shell(nullptr)
{}

GFaceUse::GFaceUse(SGModel *m, GFace * face, GLoopUse * loop, GShell *shell)
  : GEntity(m,0), Face(face), Shell(shell)
{
  Loops.append(loop);
}

GeoRep * GFaceUse::geometry()
{ 
  notImplemented();
  return nullptr;
}

TopoType::Value GFaceUse::type() const
{ 
  return TopoType::FaceUse;
}

int GFaceUse::dim() const
{ return 2; }

AOMD::SBoundingBox3d GFaceUse::bounds() const
{ return Face->bounds(); }

GShell * GFaceUse::shell() const
{ return Shell; }

GFace * GFaceUse::face() const
{ return Face; }

SSListCIter<GLoopUse*> GFaceUse::firstLoopUse() const
{ return SSListCIter<GLoopUse*>(Loops); }

int GFaceUse::dir() const
{
  return Face->whichUse(this);
}

void GFaceUse::addLoopUse(GLoopUse * lu)
  // first loopuse added must be outer loop
{ 
  Loops.append(lu); 
  lu->setFaceUse(this);
}

void GFaceUse::setShell(GShell *shell)
{ Shell = shell; }

int GFaceUse::inClosure(GEntity *ent) const
{
  if(this==ent || Face == ent)
    return 1;
  GLoopUse *lu;
  for(SSListCIter<GLoopUse*> luIter(Loops); luIter(lu); )
    if(lu->inClosure(ent))
      return 1;
  return 0;
}

SSList<GEdge *> GFaceUse::edges() const
{
  SSList<GEdge*> es;
  GLoopUse *lu;
  for(SSListCIter<GLoopUse*> luIter(Loops); luIter(lu); )
    es.appendUnique(lu->edges());
  return es;
}


GLoopUse * GFaceUse::insertEdge(GEdge *e)
  // insert a new edge into the definition of a faceuse, the new edge
  // points to the two vertices at its ends which must be on the 
  // closure of the face
  // if the inserted edge splits the face, the new loop use created for
  // the new loop of edges is returned, otherwise null is returned
{
  GLoopUse *retLoop = nullptr;
  GVertex *newv[2];
  newv[0] = e->vertex(0);
  newv[1] = e->vertex(1);
  GEdgeUse *eusBefore[2] = {nullptr,nullptr}; // edgeUseSide before vertex
  GLoopUse * loopUse[2] = {nullptr,nullptr}; // loop uses

  GLoopUse *lu;
  SSListCIter<GLoopUse*> luIter(Loops);
  while(luIter(lu)){
    GEdgeUse *eus;
    SSListCIter<GEdgeUse*> eusIter(lu->firstEdgeUse());
    while(eusIter(eus)){
      GVertex *v = eus->vertexUse()->vertex();
      for(int i = 0; i < 2; i++){
	if(v==newv[i]){
	  if(eusBefore[i] || loopUse[i])
	    cerr << "GFace::insertEdge error\n";
	  eusBefore[i] = eus;
	  loopUse[i] = lu;
	}
      }
    }
  }
  // now should have the two split points needed, check
#ifdef DEBUG
  if(eusBefore[0] == eusBefore[1])
    cerr << "GFaceUse::insertEdge - error 2\n";
  if(!eusBefore[0] || !eusBefore[1] || !loopUse[0] || !loopUse[1])
    cerr << "GFaceUse::insertEdge - didn't find insertions\n";
#endif
  if(loopUse[0] == loopUse[1]){ // splitting a loop
    // this case not done yet, want to do other one first
    int dir=1;
    retLoop = loopUse[0]->split(e, dir, eusBefore[0], eusBefore[1]);
    cerr << "GFaceUse::insertEdge - splitting face, still in testing\n";
  } else {
    // joining two loops that are already on same face
    // Note: by the way they were found, the edge is always used in
    // the positive direction for the use after eusBefore[0]
    GLoopUse *oldUse = eusBefore[1]->loopUse();
    loopUse[0]->merge(e,1,eusBefore[0],eusBefore[1]);
    Loops.remove(oldUse);
  }
  return retLoop;
}

void GFaceUse::write(ostream &out)
{
  out << Loops.size() << " ";
  GLoopUse *lu;
  for(SSListCIter<GLoopUse*> luIter(Loops); luIter(lu); )
    lu->write(out);
}

void GFaceUse::read(istream &in)
{
  int nloop;
  in >> nloop;
  for(int i = 0; i < nloop; i++){
    GLoopUse * lu = new GLoopUse(model(),this);
    lu->read(in);
    addLoopUse(lu);
  }
}


void GFaceUse::info() const
{
  cerr << "faceuse: ";
  SSListCIter<GLoopUse*> luIter(Loops);
  GLoopUse *lu;
  while(luIter(lu)){
    lu->info();
  }
  cerr << "\n";
}
