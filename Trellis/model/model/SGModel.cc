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
#include "SGModel.h"
#include "GRegion.h"
#include "GFace.h"
#include "GEdge.h"
#include "GVertex.h"
#include "GShell.h"
#include "SSList.h"
#include "TRegion.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using std::cerr;
using std::cout;


extern SGModel *theModel;

SGModel::SGModel(const SString &name)
  : SModel(name), d_Coherent(0), OuterRegion(nullptr)
{
  if(theModel)
    cerr << "Warning another model was already loaded\n";
  theModel = this;
}

SGModel::~SGModel()
{
//  int i;

//   for(i=0; i < RegionList.size(); i++)
//     delete RegionList[i];

//   for(i=0; i < FaceList.size(); i++)
//     delete FaceList[i];

//   for(i=0; i < EdgeList.size(); i++)
//     delete EdgeList[i];

//   for(i=0; i < VertexList.size(); i++)
//     delete VertexList[i];

}


SGModel * SGModel::createModelSMD(const SString &n)
{
  std::ifstream infile(n);
  SString name,type;
  int version;
  infile >> type >> version >> name;
  SGModel *model=nullptr;

  if(version > 2)
    cerr << "Unknown file version number, winging it...";

  for(int i = 0; i < Num; i++){
    if(type == Names[i])
      model = (SGModel*)(IC[i](name,version==1?0:1,nullptr));
  }
  if(!model)
    cerr << "Couldn't load model named " << name << " of type " << type << "\n";
  
  if(version==2){
    // read in rest of file
    model->readSMD(infile);
  }
  return model;
}


int SGModel::numRegion() const
{ return RegionList.size(); }

int SGModel::numFace() const
{ return FaceList.size(); }

int SGModel::numEdge() const
{ return EdgeList.size(); }

int SGModel::numVertex() const
{ return VertexList.size(); }

GRegion *SGModel::region(int n) const
{ return RegionList.nth(n); }

GFace *SGModel::face(int n) const
{ return FaceList.nth(n); }

GEdge *SGModel::edge(int n) const
{ return EdgeList.nth(n); }

GVertex *SGModel::vertex(int n) const
{ return VertexList.nth(n); }

SSListCIter<GRegion*> SGModel::firstRegion() const
{ return SSListCIter<GRegion*>(RegionList); }

SSListCIter<GFace*> SGModel::firstFace() const
{ return SSListCIter<GFace*>(FaceList); }

SSListCIter<GEdge*> SGModel::firstEdge() const
{ return SSListCIter<GEdge*>(EdgeList); }

SSListCIter<GVertex*> SGModel::firstVertex() const
{ return SSListCIter<GVertex*>(VertexList); }

GRegion * SGModel::outerRegion() const
{ return OuterRegion; }

void SGModel::stats()
{

  cout << "\tNumber of Regions: " << RegionList.size() << "\n";
  cout << "\tNumber of Faces: " << FaceList.size() << "\n";
  cout << "\tNumber of Edges: " << EdgeList.size() << "\n";
  cout << "\tNumber of Verticies: " << VertexList.size() << "\n";

}

int abssign(int *);

void SGModel::writeSMD(const SString &n)
{
  int i;
  AOMD::SBoundingBox3d bbox;
  
  std::ofstream outfile(n);
  
  outfile << modeler() << " 3\n";
  outfile << name() << "\n";
  
  // Write out the bounds of the model
  bbox = bounds();
  outfile<< bbox.min().x()<<" "<< bbox.min().y()<<" "<< bbox.min().z() <<"\n"; 
  outfile<< bbox.max().x()<<" "<< bbox.max().y()<<" "<< bbox.max().z() <<"\n";
 
  outfile << numRegion() << " " << numFace()<<" "<<numEdge()<<" "<<numVertex()<<"\n";

  i=0;
  SSListCIter<GVertex*> vIter = firstVertex();
  GVertex *v;
  while( vIter(v) ){
    v->write(outfile);
    v->setID(i++);
  }

  i=0;
  SSListCIter<GEdge*> eIter = firstEdge();
  GEdge *e;
  while( eIter(e) ){
    e->write(outfile);
    e->setID(i++);
  }

  i=0;
  SSListCIter<GFace*> fIter = firstFace();
  GFace *f;
  while( fIter(f) ){
    f->write(outfile);
    f->setID(i++);
  }

  i=0;
  SSListCIter<GRegion*> rIter = firstRegion();
  GRegion *r;
  while( rIter(r) ){
    r->write(outfile);
    r->setID(i++);
  }
  OuterRegion->write(outfile);

  outfile.close();
  }

void SGModel::readBounds(std::istream &in) {
  // Reads the min and max pnts of the model bounds from the smd file
  // and then ignores it since most of the derivied classes ignore it

  AOMD::SPoint3 pnt;
  in >> pnt;
  in >> pnt;
}

void SGModel::readSMD(std::istream &in)
{
  int i;
  SString temp1;
  GVertex *vtx;
  GEdge *ed;
  GFace *fc;
  GRegion *reg;
  
  int nr, nf, ne, nv;
  
  // Skip modeler name and read in version
  in >> temp1 >> i;
  
  if (i != 3)
    cerr << "SMD File is a Version " << i 
	 << " File that is not supported - this may cause problems\n";

  // Read in model name
  in >> Name;

  // Read in bounds
  readBounds(in);

  in >> nr >> nf >> ne >> nv;

  for(i=0; i < nv; i++){
    vtx = createVertex(in);
    add(vtx);
    vtx->read(in);
  }

  for(i=0; i < ne; i++){
    ed = createEdge(in);
    add(ed);
    ed->read(in);
  }

  for(i=0; i < nf; i++){
    fc = createFace(in);
    add(fc);
    fc->read(in);
  }

  for(i=0; i < nr; i++){
    reg = createRegion(in);
    add(reg);
    reg->read(in);
  }
  int tag;
  in >> tag;
  OuterRegion = new TRegion(this,tag);
  OuterRegion->read(in);
}


GRegion * SGModel::regionByTag(int tag) const
{
  GRegion * e;
  SSListCIter<GRegion*> eIter = firstRegion();
  while(eIter(e)){
    if(e->tag() == tag)
      return e;
  }
  return nullptr;
}

GFace * SGModel::faceByTag(int tag) const
{
  GFace * e;
  SSListCIter<GFace*> eIter = firstFace();
  while(eIter(e)){
    if(e->tag() == tag)
      return e;
  }
  return nullptr;
}

GEdge * SGModel::edgeByTag(int tag) const
{
  GEdge * e;
  SSListCIter<GEdge*> eIter = firstEdge();
  while(eIter(e)){
    if(e->tag() == tag)
      return e;
  }
  return nullptr;
}


GVertex * SGModel::vertexByTag(int tag) const
{
  GVertex * e;
  SSListCIter<GVertex*> eIter = firstVertex();
  while(eIter(e)){
    if(e->tag() == tag)
      return e;
  }
  return nullptr;
}


GRegion * SGModel::regionByID(int ident) const
{
  SSListCIter<GRegion*> rIter = firstRegion();
  GRegion *r;
  while( rIter(r) ){
    if( r->id() == ident)
      return r;
  }
  return nullptr;
}

GFace * SGModel::faceByID(int ident) const
{
  SSListCIter<GFace*> fIter = firstFace();
  GFace *f;
  while( fIter(f) ){
    if( f->id() == ident )
      return f;
  }
  return nullptr;

}

GEdge * SGModel::edgeByID(int ident) const
{
  SSListCIter<GEdge*> eIter = firstEdge();
  GEdge *e;
  while( eIter(e) ){
    if( e->id() == ident )
      return e;
  }
  return nullptr;

}

GVertex * SGModel::vertexByID(int ident) const
{
  SSListCIter<GVertex*> vIter = firstVertex();
  GVertex *v;
  while( vIter(v) ){
    if( v->id() == ident )
      return v;
  }
  return nullptr;

}

SModelMember * SGModel::getMember(int type, int tag) const
{
  switch(type){
  case TopoType::Vertex:
    return vertexByTag(tag);
  case TopoType::Edge:
    return edgeByTag(tag);
  case TopoType::Face:
    return faceByTag(tag);
  case TopoType::Region:
    return regionByTag(tag);
  case TopoType::Model:
    return (SModelMember*)this; // cast away const
  default:
    cerr << "SGModel::getMember - error, type = " << type << "\n";
    break;
  }
  return nullptr;
}

SModelMember * SGModel::getMember(const SString &name) const
{
  // should really use a hash table for this
  int i;
  for(i=0; i < numRegion(); i++)
    if( region(i)->name() == name )
      return region(i);

  for(i=0; i < numFace(); i++)
    if( face(i)->name() == name )
      return face(i);

  for(i=0; i < numEdge(); i++)
    if( edge(i)->name() == name )
      return edge(i);

  for(i=0; i < numVertex(); i++)
    if( vertex(i)->name() == name )
      return vertex(i);
  return nullptr;

}

void SGModel::add(GRegion *r)
{ 
  r->setID(numRegion());
  RegionList.append(r); 
}

void SGModel::add(GFace *f)
{ 
  f->setID(numFace());
  FaceList.append(f); 
}

void SGModel::add(GEdge *e)
{ 
  e->setID(numEdge());
  EdgeList.append(e); 
}

void SGModel::add(GVertex *v)
{ 
  v->setID(numVertex());
  VertexList.append(v); 
}

void SGModel::remove(GRegion *r)
{ RegionList.remove(r); }

void SGModel::remove(GFace *f)
{ FaceList.remove(f); }

void SGModel::remove(GEdge *e)
{ EdgeList.remove(e); }

void SGModel::remove(GVertex *v)
{ VertexList.remove(v); }

void SGModel::createOuterShell()
{
  SSList<GFace*> faces;
  SSList<int> dirs;
  for(int i=0; i < numFace(); i++){
    GFace *f = face(i);
    for(int dir =0; dir < 2; dir++){
      if(f->use(dir)->shell() == nullptr){
	faces.append(f);
	dirs.append(dir);
      }
    }
  }
  OuterRegion = new TRegion(this,0,nullptr,faces,dirs);
  //new GShell(r,faces,dirs); // should do something with this
}

void SGModel::setDisplayCoherence(int c)
{ d_Coherent = c; }


