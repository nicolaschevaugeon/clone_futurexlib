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
#include "GRegion.h"
#include "GShell.h"
#include "GFaceUse.h"
#include "GFace.h"
#include "GEdge.h"
#include "GVertex.h"
#include "SSList.h"
#include "SGModel.h"
#include "MessageOut.h"

GRegion::GRegion(SGModel *m, int t, const SSList<GFace*> & faces, const SSList<int> &dirs)
  : GEntity(m,t)
{
  GShell *su = new GShell(m,this,faces,dirs);
  Shells.append(su);
}

GRegion::GRegion(SGModel *m, int t)
  : GEntity(m, t)
{
  //cerr << "GRegion constructor 2 called for " << this << endl;
}

GRegion::~GRegion() 
{
  // add destruction of all lists should happen by SSList destructors.
}

TopoType::Value GRegion::type() const
{ return TopoType::Region; }

int GRegion::dim() const
{ return 3; }

SSList<GRegion*> GRegion::regions() const
{
  SSList<GRegion*> rs;
  return rs;
}

SSList<GFace*> GRegion::faces() const
{
  SSList<GFace*> fs;
  GShell *su;
  for(SSListCIter<GShell*> suIter = firstShell(); suIter(su); ){
    GFaceUse *fu;
    for(SSListCIter<GFaceUse*> fuIter = su->firstFaceUse(); fuIter(fu); )
      fs.appendUnique(fu->face());
  }
  return fs;
}

SSList<GEdge*> GRegion::edges() const
{
  SSList<GEdge*> es;
  GShell *su;
  for(SSListCIter<GShell*> suIter = firstShell(); suIter(su); ){
    GFaceUse *fu;
    for(SSListCIter<GFaceUse*> fuIter = su->firstFaceUse(); fuIter(fu); )
      es.appendUnique(fu->face()->edges());
  }

  return es;
}

SSList<GVertex*> GRegion::vertices() const
{
  SSList<GVertex*> vs;
  GShell *su;
  for(SSListCIter<GShell*> suIter = firstShell(); suIter(su); ){
    GFaceUse *fu;
    for(SSListCIter<GFaceUse*> fuIter = su->firstFaceUse(); fuIter(fu); )
      vs.appendUnique(fu->face()->vertices());
  }

  return vs;
}

int GRegion::classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const
{
  if(chkb){
    GVertex *v;
    for (SSListCIter<GVertex *> vIter(vertices()); vIter(v); ){
      if(v->containsPoint(pt)){
	if(bent)
	  *bent = v;
	return -1;
      }
    }
    GEdge *e;
    for (SSListCIter<GEdge *> eIter(edges()); eIter(e); ){
      if(e->containsPoint(pt)){
	if(bent)
	  *bent = e;
	return -1;
      }
    }
    GFace *f;
    for (SSListCIter<GFace *> fIter(faces()); fIter(f); ){
      if(f->containsPoint(pt)){
	if(bent)
	  *bent = f;
	return -1;
      }
    }
  }
  if(containsPoint(pt))
    return 1;
  return 0;
}

SSListCIter<GShell*> GRegion::firstShell() const
{
  return SSListCIter<GShell*>(Shells);
}

AOMD::SBoundingBox3d GRegion::bounds() const
// for a region, same as first (outer) shell
{ 
  return Shells.nth(0)->bounds();
}

void GRegion::appendShell(GShell *shell)
{ 
  Shells.append(shell); 
  shell->setRegion(this);
}

void GRegion::addShell(const SSList<GFace*> & faces, const SSList<int> &dirs)
{
  GShell *su = new GShell(model(),this,faces,dirs);
  Shells.append(su);
}

int GRegion::inClosure(GEntity *ent) const
{
  if(this==ent)
    return 1;
  GShell *su;
  for(SSListCIter<GShell*> suIter(Shells); suIter(su); )
    if(su->inClosure(ent))
      return 1;
  return 0;
}

void GRegion::write(ostream &out)
{
  out << tag() << " "  << id() << " " << Shells.size() << " \n";
  GShell *s;
  for(SSListCIter<GShell *> suIter(Shells); suIter(s); )
    s->write(out);
  out << "\n";
}

void GRegion::read(istream &in)
{
  int nshell, id;
  in >> id >> nshell;
  setID(id);
  for(int i =0; i < nshell; i++){
    GShell *s = new GShell(model(),this);
    s->read(in);
    appendShell(s);
  }
}




