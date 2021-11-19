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
#include "NullModel.h"

int NullModel::nullSetup = 0;

#include <iostream>
#include <cstdlib>

#include "NullVertex.h"
#include "NullEdge.h"
#include "NullFace.h"
#include "NullRegion.h"
#include "GLoopUse.h"
#include "GShell.h"
#include "SString.h"
#include "SBoundingBox3d.h"

extern "C"
void GM_registerNull()
{ NullModel::registerModeler(); }

extern "C"
SGModel * GM_createFromNothing(const char *name)
{ return new NullModel(name); }

NullModel::NullModel()
: SGModel("unnamed")
{
  NullModel::setupNull();
}

void NullModel::registerModeler()
{ SGModel::registerModeler("null",NullModel::create); }

SModel * NullModel::create(const SString &name, int loadOnly, 
			    const SSList<SModel*> *const)
{ return new NullModel(name); }


NullModel::NullModel(const SString &name)
  : SGModel(name)
{
  NullModel::setupNull();

}

SString NullModel::modeler() const
{ return SString("null"); }

AOMD::SBoundingBox3d NullModel::bounds() const
{
  AOMD::SBoundingBox3d box;
  return box;
}

double NullModel::tolerance() const
{
  return 1e-6;;
}

void NullModel::setGeomTolerance(double tol)
{
}


GVertex * NullModel::createVertex(istream &in)
{
  int tag;
  in >> tag;
  return new NullVertex(this,tag);
}

GEdge * NullModel::createEdge(istream &in)
{
  int tag;
  in >> tag;
  return new NullEdge(this,tag);
}

GFace * NullModel::createFace(istream &in)
{
  int tag;
  in >> tag;
  return new NullFace(this,tag);
}

GRegion * NullModel::createRegion(istream &in)
{
  int tag;
  in >> tag;
  return new NullRegion(this,tag);
}

void NullModel::buildFromGeom()
{
}

extern "C" 
void GM_setupNull()
{
  NullModel::setupNull();
}

void NullModel::setupNull()
{
  if(!NullModel::nullSetup){
    NullModel::nullSetup = 1;
  }
}

GRegion * NullModel::regionByTag(int tag) const
{
  GRegion *e = SGModel::regionByTag(tag);
  if(!e){
    e = new NullRegion((NullModel*)this,tag);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;
}

GFace * NullModel::faceByTag(int tag) const
{
  GFace *e = SGModel::faceByTag(tag);
  if(!e){
    e = new NullFace((NullModel*)this,tag);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;

}

GEdge * NullModel::edgeByTag(int tag) const
{
  GEdge *e = SGModel::edgeByTag(tag);
  if(!e){
    e = new NullEdge((NullModel*)this,tag);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;

}

GVertex * NullModel::vertexByTag(int tag) const
{
  GVertex *e = SGModel::vertexByTag(tag);
  if(!e){
    e = new NullVertex((NullModel*)this,tag);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;

}

GRegion * NullModel::regionByID(int ident) const
{
  GRegion *e = SGModel::regionByID(ident);
  if(!e){
    e = new NullRegion((NullModel*)this,ident);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;

}

GFace * NullModel::faceByID(int ident) const
{
  GFace *e = SGModel::faceByID(ident);
  if(!e){
    e = new NullFace((NullModel*)this,ident);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;

}

GEdge * NullModel::edgeByID(int ident) const
{
  GEdge *e = SGModel::edgeByID(ident);
  if(!e){
    e = new NullEdge((NullModel*)this,ident);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;

}

GVertex * NullModel::vertexByID(int ident) const
{
  GVertex *e = SGModel::vertexByID(ident);
  if(!e){
    e = new NullVertex((NullModel*)this,ident);
    ((NullModel*)this)->add(e); // cast away const
  }
  return e;

}
