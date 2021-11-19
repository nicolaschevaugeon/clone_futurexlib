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
#include "GShell.h"
#include "GRegion.h"
#include "GFace.h"
#include "GEdge.h"
#include "GVertex.h"
#include "SGModel.h"
#include "SSList.h"
using std::cerr;

GShell::GShell(SGModel *m,GRegion *region, const SSList<GFace*> &faces, const SSList<int> &dirs)
  : GEntity(m,0), Region(region)
{
  int dir;
  GFace *face;
  SSList<GEdge*> es;
  SSListCIter<GFace *> fIter(faces);
  SSListCIter<int> dIter(dirs);
  for( ; fIter(face), dIter(dir) ; ){
    GFaceUse *fu = appendFace(face,dir);
    es.appendUnique(fu->edges());
  }
  GEdge *e;
  for(SSListIter<GEdge*> eIter(es); eIter(e); )
    e->mergeUses(this);
}

GShell::GShell(SGModel *m,GRegion *region)
  : GEntity(m,0), Region(region)
{
}

GShell::~GShell()
= default;

GeoRep * GShell::geometry()
{ 
  cerr << "GShell::geometry() called\n"; 
  return nullptr;
}

//void GShell::setMesh(MEntity *ent)
//{ cerr << "GShell::setMesh() called\n"; }

TopoType::Value GShell::type() const
{ return TopoType::Shell; }

int GShell::dim() const
{ return 2; }

GRegion* GShell::region() const
{ 
  return (Region == model()->outerRegion() ? nullptr : Region); 
}

SSListCIter<GFaceUse*> GShell::firstFaceUse() const
{ return SSListCIter<GFaceUse*>(FaceUses);}

AOMD::SBoundingBox3d GShell::bounds() const
// for a shell bounding box is bounds of all faces
{
  AOMD::SBoundingBox3d box;
  GFaceUse *f;
  for(SSListCIter<GFaceUse*> fIter(FaceUses); fIter(f); )
    box += f->bounds();
  return box;
}

void GShell::setRegion(GRegion *region)
{ Region = region; }

GFaceUse * GShell::appendFace(GFace *face, int dir)
{
  GFaceUse *fu = face->use(dir);
  FaceUses.append(fu);
  fu->setShell(this);
  return fu;
}

int GShell::inClosure(GEntity *ent) const
{
  if(this==ent)
    return 1;
  GFaceUse *fu;
  for(SSListCIter<GFaceUse*> fuIter(FaceUses); fuIter(fu); )
    if(fu->inClosure(ent))
      return 1;
  return 0;
}

void GShell::removeFaceUse(GFaceUse *fu)
{  FaceUses.remove(fu); }

void GShell::write(ostream &out)
{
  GFaceUse *fu;
  out << FaceUses.size() << "\n";
  for(SSListCIter<GFaceUse*> fuIter(FaceUses); fuIter(fu); ){
    out << fu->face()->id() << " " << fu->dir() << " ";
  }
}

void GShell::read(istream &in)
{
  int nface;
  in >> nface;
  for(int i = 0; i < nface; i++){
    int faceid, dir;
    in >> faceid >> dir;
    GFaceUse *fu = model()->faceByID(faceid)->use(dir);
    FaceUses.append(fu);
    fu->setShell(this);
  }
}
