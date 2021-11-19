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
#include "TFace.h"
#include "SVector3.h"
#include <iostream>
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
#include "GRegion.h"
#include "GFPoint.h"
#include "GEPoint.h"
#include "GVPoint.h"
#include "TEdge.h"
#include "GEdgeUsePair.h"
#include "GVertexUse.h"
#include "TFaceSplitEdge.h"
#include "TPartialFace.h"
#include "GShell.h"
#include "SGModel.h"
#include "SSList.h"

using std::cerr;

TFace::TFace(SGModel *model, int t, GFace *f,SSList<GEdge*> edges, 
	     SSList<int> dirs)
  : GFace(model, t, edges,dirs), RealFace(f)
{
  if(f){
    d_periodic[0] = f->periodic(0);
    d_periodic[1] = f->periodic(1);
  }
}

TFace::~TFace()
= default;

AOMD::Range<double> TFace::parBounds(int i) const
{ 
  if (RealFace)
    return RealFace->parBounds(i); 
  return AOMD::Range<double>(-1.0e100, 1.0e100);
}

int TFace::paramDegeneracies(int dir, double *par)
{ return RealFace->paramDegeneracies(dir,par); }

AOMD::SBoundingBox3d TFace::bounds() const
{ return RealFace->bounds(); }

SVector3 TFace::normal(const SPoint2 &param) const
{ return RealFace->normal(param);}

Pair<SVector3,SVector3> TFace::firstDer(const SPoint2 &param) const
{ return RealFace->firstDer(param); }

double * TFace::nthDerivative(const SPoint2 &param, int n, double *array) const
{ return RealFace->nthDerivative(param,n,array); }

GFPoint TFace::point(const SPoint2 &pt) const
{ return RealFace->point(pt); }

GFPoint TFace::point(double par1,double par2) const
{ return RealFace->point(par1,par2); }

GFPoint TFace::closestPoint(const AOMD::SPoint3 & queryPoint)
{ return RealFace->closestPoint(queryPoint); }

int TFace::containsPoint(const AOMD::SPoint3 &pt) const
{ return RealFace->containsPoint(pt); }

int TFace::containsParam(const SPoint2 &par) const
{ return RealFace->containsParam(par); }

SPoint2 TFace::parFromPoint(const AOMD::SPoint3 &pt) const
{ return RealFace->parFromPoint(pt); }

GeoRep * TFace::geometry()
{ return RealFace->geometry(); }

Logical::Value TFace::continuous(int dim) const
{ return RealFace->continuous(dim); }

Logical::Value TFace::periodic(int dim) const
{ 
  if(d_periodic[dim])
    return Logical::True;
  else
    return Logical::False;
}

Logical::Value TFace::degenerate(int dim) const
{ return RealFace->degenerate(dim); }

GeomType::Value TFace::geomType() const
{ return RealFace->geomType(); }

int TFace::geomDirection() const
{ return RealFace->geomDirection(); }

double TFace::tolerance() const
{ 
  if (RealFace)
    return RealFace->tolerance(); 
  return 1.0e-6;
}


SSList<GEdge*> TFace::splitU(double u)
{
  return splitAtConstParam(0,u,u);
}

SSList<GEdge*> TFace::splitV(double v)
{ return splitAtConstParam(1,v,v); }

SSList<GEdge*> TFace::splitSeamU(double ulow, double uhigh)
{
  return splitAtConstParam(0,ulow,uhigh);
}

SSList<GEdge*> TFace::splitSeamV(double vlow, double vhigh)
{ return splitAtConstParam(1,vlow,vhigh); }

SSList<GEdge*> TFace::splitAtConstParam(int axis, double param1, double param2)
  // axis is the parametric direction being split
  // param1 and param2 are the values in the axis direction along
  // each side of the edge.

  // Note: this really can't be a member function of the face, since
  // the face being split could disappear before the function finishes
  // should probably be a member function of the model.
{
  // in some cases this might not actually create a new face (cylinder
  // with two "outer" loops

  // first need to get a list of vertices along the split line
  // if they are already there, use them, if not, make them
  GVertex *intsV[10];
  double intsPar[10];
  int nint = 0;
  
  GEdge *edge;
  for(SSListCIter<GEdge*> eIter(edges()); eIter(edge); ){
    SSList<GEPoint> epts = edge->intersect(axis,param1,this);
    GEPoint ept;
    // loop over all intersection points and make a vertex there
    // if needed by spliting the edge

    // NOTE: When there are multiple splits of an edge, some of the split
    // points may be on newly generated edges after one or more of 
    // the splits is done, this is not handled correctly here
    for(SSListIter<GEPoint> eptIter(epts); eptIter(ept); ){
      GVertex *ev= edge->vertex(0)->containsPoint(ept) ? edge->vertex(0) : nullptr;
      if(!ev)
	ev = edge->vertex(1)->containsPoint(ept) ? edge->vertex(1) : nullptr;
      if(ev) { // have vertex
	for(int i=0; i < nint; i++){
	  if(intsV[i] == ev){
	    ev = nullptr;
	    break;
	  }
	}
	if(ev){
	  intsV[nint] = ev;
	  SPoint2 fpt = edge->reparamOnFace(this,ept.par(),1); //need right direction?
	  intsPar[nint] = axis ? fpt[0] : fpt[1];
	  nint++;
	}
      } else {
	cerr << "split edge\n";
	SPoint2 fpt = edge->reparamOnFace(this,ept.par(),1); // direction?
	intsV[nint] = edge->split(ept.par());
	intsPar[nint] = axis ? fpt[0] : fpt[1];
	nint++;
      }
    }
  }
  cerr << "found " << nint << " splits\n";

  SSList<GEdge*> newEdges;

  // sort in order of increasing parameter location, badly written 
  // bubble sort
  int done = 0;
  while(!done){
    done = 1;
    for(int i = 0; i < nint-1; i++){
      if(intsPar[i] > intsPar[i+1]){
	done = 0;
	double tpar = intsPar[i+1];
	GVertex *tv = intsV[i+1];
	intsPar[i+1] = intsPar[i];
	intsV[i+1] = intsV[i];
	intsPar[i] = tpar;
	intsV[i] = tv;
      }
    }
  }

  // first just deal with case of a single new edge
  if(nint != 2)
    cerr << "not two intersections\n";

  GEdge *e;
  //int edir;
  e = new TFaceSplitEdge(model(),RealFace->tag(),RealFace,intsV[0],intsV[1],
			 intsPar[0],intsPar[1],axis,param2,param1);
  newEdges.append(e);
  model()->add(e);

  if(d_periodic[axis])
    d_periodic[axis]=0;

  GLoopUse * luf = FrontUse.insertEdge(e);
  GLoopUse * lub = BackUse.insertEdge(e);  
  if(luf && lub){
    //TPartialFace(SGModel *model, GFace *f, int axis, int side, double par);
    TPartialFace *newFace = new TPartialFace(model(),RealFace->tag(),RealFace,
					     axis,1,param2,
					     luf,FrontUse.shell(),
					     lub,BackUse.shell());
  }


  return newEdges;

}

void TFace::makeNewLoops(GEdge *e, GFaceUse *fu, GFaceUse *fu1, GFaceUse *fu2, 
			 GLoopUse **l1, GLoopUse **l2)
  // given fu split it into two by introducing e
{
  // e is the new edge
  GVertex *v0 = e->vertex(0);
  GVertex *v1 = e->vertex(1);
  SSList<GEdgeUse*> eu1;
  SSList<GEdgeUse*> eu2;
  GEdgeUsePair *newUse;
  int loopNum = 0;
  SSListCIter<GLoopUse*> luIter = fu->firstLoopUse();
  GLoopUse *lu;
  luIter(lu);
  GEdgeUse *eu;
  for(SSListCIter<GEdgeUse*> euIter = lu->firstEdgeUse(); euIter(eu); ){
    int euDir = eu->dir();
    GVertexUse *vu= eu->vertexUse(); // VU at head of EU
    if(loopNum==0 && vu->vertex() == v0){
      loopNum = 2;
      newUse = e->addLoop(nullptr,1)->use();
      // add use of e to eu1 (+ve dir)
      // add use of e to eu2 (-ve dir)
      continue;
    }
    if(loopNum == 2 && vu->vertex() == v1){
      eu2.append(eu);  // done with loop in eu2
      loopNum = 1;
      continue;
    }
    if(loopNum==1)
      eu1.append(eu);
    else if (loopNum ==2)
      eu2.append(eu);
  }
  newUse->side(1)->setLoopUse(lu);
  *l1 = new GLoopUse(model(), eu1,fu1,lu);  // this replaces old edgeUse-loopUse ptrs
  newUse->side(0)->setLoopUse(lu);
  *l2 = new GLoopUse(model(),eu2,fu2,lu);
}

AOMD::GEntity * TFace::getBaseRep()
{ return RealFace->getBaseRep(); }

void * TFace::getNativePtr() const
{ return RealFace->getNativePtr(); }

int TFace::getNativeInt() const
{ return RealFace->getNativeInt(); }

Logical::Value TFace::surfPeriodic(int dim) const
{ return RealFace->surfPeriodic(dim); }

double TFace::period(int dir) const
{ return RealFace->period(dir); }

void TFace::insertEdge(GEdge *e)
  // insert a new edge into the definition of a face, the new edge
  // points to the two vertices at its ends which must be on the 
  // closure of the face
  // this may split the face, returns a newly created face if this happens
{

}

//  Attribute * TFace::attribute(const SString &type) const
//  {
//    if(RealFace)
//      return RealFace->attribute(type);
//    return  SModelMember::attribute(type); 
//  }

//  SSList<Attribute *> TFace::attributes(const SString &type) const
//  {
//    if(RealFace)
//      return RealFace->attributes(type);
//    return  SModelMember::attributes(type); 
//  }

//  SSList<Attribute *> TFace::attributes() const
//  {
//    if(RealFace)
//      return RealFace->attributes();
//    return  SModelMember::attributes(); 
//  }
