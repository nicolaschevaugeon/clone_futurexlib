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
#include "TopoModel.h"
#include "TVertex.h"
#include "TEdge.h"
#include "TFace.h"
#include "TRegion.h"
#include "GEdgeUsePair.h"
#include "GLoopUse.h"
#include "GFaceUse.h"
#include "GShell.h"

using std::cerr;

// create a topological copy of an existing model
TopoModel::TopoModel(SGModel *source)
:  SGModel(source->name()), RealModel(source)
{

  GVertex *v;
  for(SSListCIter<GVertex*> vIter(source->firstVertex()); vIter(v); ){
    TVertex *tv = new TVertex(this,v->tag(),v);
    v->setDataP("newe",tv);
    add(tv);
  }

  GEdge *e;
  for(SSListCIter<GEdge*> eIter(source->firstEdge()); eIter(e); ){
    TEdge *te = new TEdge(this,e->tag(),e,(GVertex*)e->vertex(0)->dataP("newe"),
			  (GVertex*)e->vertex(1)->dataP("newe"));
    e->setDataP("newe",te);
    add(te);
  }

  GFace *f;
  for(SSListCIter<GFace*> fIter(source->firstFace()); fIter(f); ){
    SSList<GEdge *> edges;
    SSList<int> dirs;
    GFaceUse *fu = f->use(1);
    SSListCIter<GLoopUse*> luIter = fu->firstLoopUse();
    GLoopUse *lu;
    luIter(lu); // get outer loop
    GEdgeUse *eu;
    for(SSListCIter<GEdgeUse*> euIter = lu->firstEdgeUse(); euIter(eu); ){
      edges.append((GEdge*)eu->edge()->dataP("newe"));
      dirs.append(eu->dir());
    }
    TFace *tf = new TFace(this,f->tag(),f,edges,dirs);
    f->setDataP("newe",tf);
    while(luIter(lu)){  // inner loops
      edges.clear();
      dirs.clear();
      GEdgeUse *eui;
      for(SSListCIter<GEdgeUse*> euiIter = lu->firstEdgeUse(); euiIter(eui); ){
	edges.append((GEdge*)eui->edge()->dataP("newe"));
	dirs.append(eui->dir());
      }
      tf->addLoop(edges,dirs);
    }
    add(tf);
  }

  GRegion *r;
  for(SSListCIter<GRegion*> rIter(source->firstRegion()); rIter(r); ){
    SSList<GFace *> faces;
    SSList<int> dirs;
    SSListCIter<GShell*> suIter = r->firstShell();
    GShell *su;
    suIter(su); // get outer shell
    GFaceUse *fu;
    for(SSListCIter<GFaceUse*> fuIter = su->firstFaceUse(); fuIter(fu); ){
      faces.append((GFace*)fu->face()->dataP("newe"));
      dirs.append(fu->dir());
    }
    TRegion *tf = new TRegion(this,r->tag(),r,faces,dirs);
    while(suIter(su)){  // inner shells
      faces.clear();
      dirs.clear();
      GFaceUse *fui;
      for(SSListCIter<GFaceUse*> fuiIter = su->firstFaceUse(); fuiIter(fui); ){
	faces.append((GFace*)fui->face()->dataP("newe"));
	dirs.append(fui->dir());
      }
      cerr << "TopoModel::TopoModel - inner shell on region - not handled yet\n";
    }
    add(tf);
  }
  createOuterShell();
}

TopoModel::TopoModel(const SString &name, SGModel *source)
:  SGModel(name), RealModel(source)
{
}


SString TopoModel::modeler() const
{ return "topo"; }

AOMD::SBoundingBox3d TopoModel::bounds() const
{ 
  if(RealModel)
    return RealModel->bounds(); 

  return bbox;
    
}

void TopoModel::setBounds(const AOMD::SBoundingBox3d &box)
{ 
  bbox = box;
}

double TopoModel::tolerance() const
{
  if (RealModel)
    return RealModel->tolerance(); 
  return 1.0e-6;
}

GVertex * TopoModel::createVertex(istream &in)
{
  int tag;
  in >> tag;
  return new TVertex(this, tag, nullptr);
}

GEdge * TopoModel::createEdge(istream &in)
{
  int tag;
  in >> tag;
  return new TEdge(this, tag, nullptr, nullptr, nullptr);
}

GFace * TopoModel::createFace(istream &in)
{
  int tag;
  SSList<GEdge *> edges;
  SSList<int> dirs;
  in >> tag;
  return new TFace(this, tag, nullptr, edges, dirs);
}

GRegion * TopoModel::createRegion(istream &in)
{
  int tag;
  SSList<GFace *> faces;
  SSList<int> dirs;
  in >> tag;
  return new TRegion(this, tag, nullptr, faces, dirs);
}

void TopoModel::finalize()
{
  createOuterShell();
}

void TopoModel::augmentPeriodicFaces()
{
  // have to repeated loop over all the face until none need to
  // be split since the spliting of a face modifies the list
  // being iterated over. This isn't very optimal and should be
  // fixed
  GFace *f;
  int done = 0;
  while(!done){
    done = 1;
    for(SSListCIter<GFace*> fIter(firstFace()); fIter(f); ){
      if(!f->periodic(0) && !f->periodic(1))
	continue;
      done = 0;
      int splitParam = f->periodic(0) ? 0 : 1;
      
      // assume that the parametric jump is the extremes of the param. bounds
      AOMD::Range<double> splitRange = f->parBounds(splitParam);
      
      if(splitParam==0)
	f->splitSeamU(splitRange.low(),splitRange.high());
      else
	f->splitSeamV(splitRange.low(),splitRange.high());
      break;
    }
  }
}

//  Attribute * TopoModel::attribute(const SString &type) const
//  {
//    if(RealModel)
//      return RealModel->attribute(type);
//    return  SModelMember::attribute(type); 
//  }

//  SSList<Attribute *> TopoModel::attributes(const SString &type) const
//  {
//    if(RealModel)
//      return RealModel->attributes(type);
//    return  SModelMember::attributes(type); 
//  }

//  SSList<Attribute *> TopoModel::attributes() const
//  {
//    if(RealModel)
//      return RealModel->attributes();
//    return  SModelMember::attributes(); 
//  } 

// Read in the bounds of a model in a smd file
void TopoModel::readBounds(istream &in) {
  AOMD::SPoint3 pnt;
  in >> pnt;
  bbox += pnt;
  in >> pnt;
  bbox += pnt;
}
