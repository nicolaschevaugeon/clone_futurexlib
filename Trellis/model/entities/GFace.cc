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
#include "GFace.h"
#include "GShell.h"
#include "GFPoint.h"
#include "GEPoint.h"
#include "GVPoint.h"
#include "GLoopUse.h"
#include "GEdgeUsePair.h"
#include "GEdge.h"
#include "GVertex.h"
#include "GVertexUse.h"
#include "SSList.h"
#include "SGModel.h"
#include "modeler.h"
#include "SString.h"

using std::cout;
using std::cerr;


GFace::GFace(SGModel *m, int t, SSList<GEdge*> edges, SSList<int> dirs)
: GEntity(m,t), FrontUse(m,this), BackUse(m,this)
{
  if(edges.size()){
    // *** this should just be a call to addLoop() ***
    // need to make two loop uses
    GLoopUse *l1 = new GLoopUse(m, &FrontUse,edges,dirs);
    FrontUse.addLoopUse(l1);
    edges.reverse();
    dirs.reverse();
    int d;
    for(SSListIter<int> dIter(dirs); dIter(d); )
      dIter.replaceCurrent(d?0:1);
    GLoopUse *l2 = new GLoopUse(m,&BackUse,edges,dirs);
    BackUse.addLoopUse(l2);
    l1->setMate(l2);
    l2->setMate(l1);
  }
}

GFace::GFace(SGModel *m, int t)
  : GEntity(m,t), FrontUse(m,this), BackUse(m,this)
{ 
}

GFace::GFace(SGModel *model, int t, GLoopUse *frontLoop, GShell *frontShell,
	     GLoopUse *backLoop, GShell *backShell)
  : GEntity(model, t), FrontUse(model,this,frontLoop,frontShell),
    BackUse(model,this,backLoop,backShell)
{ 
}

GFace::~GFace()
{ 

  //delete MeshFaces;
  // need to delete the actual lists
}

TopoType::Value GFace::type() const
{ return TopoType::Face; }

int GFace::dim() const
{ return 2; }

SSList<GRegion*> GFace::regions() const
{
  SSList<GRegion*> rs;
  GRegion *r;
  if((r=region(0)))
    rs.append(r);
  if((r=region(1)))
    rs.append(r);
  return rs;
}

SSList<GFace*> GFace::faces() const
{
  SSList<GFace*> fs;
  return fs;
}

SSList<GEdge*> GFace::edges() const
{
  SSList<GEdge*> es;
  const GFaceUse *fu=use(1);
  GLoopUse *lu;
  for(SSListCIter<GLoopUse*> luIter = fu->firstLoopUse(); luIter(lu); ){
    GEdgeUse *eu;
    for(SSListCIter<GEdgeUse*> euIter = lu->firstEdgeUse(); euIter(eu); )
      es.appendUnique(eu->edge());
  }

  return es;
}

SSList<GVertex*> GFace::vertices() const
{
  SSList<GVertex*> vs;
  const GFaceUse *fu=use(1);
  GLoopUse *lu;
  for(SSListCIter<GLoopUse*> luIter = fu->firstLoopUse(); luIter(lu); ){
    GEdgeUse *eu;
    for(SSListCIter<GEdgeUse*> euIter = lu->firstEdgeUse(); euIter(eu); )
      vs.appendUnique(eu->edge()->vertices());
  }
  return vs;
}


const GFaceUse * GFace::use(int dir) const
{  return dir ? &FrontUse : &BackUse; }

GFaceUse * GFace::use(int dir)
{ return dir ? &FrontUse : &BackUse; }

GRegion * GFace::region(int dir) const
{
  GRegion *r = nullptr;
  GShell *s =  use(dir)->shell();
  if(s)
    r = s->region();
  return r;
}

SPoint2 GFace::param(const GPoint & pt) const
{
  switch( pt.gEnt()->dim() ){
  case 0:
    return parFromPoint(pt);
  case 1:{
    GEdge *e = (GEdge*)(pt.gEnt());
    return e->reparamOnFace( (GFace*)this, ((GEPoint&)pt).par(),0 );
  }
  default:
    cout << "Face::par() unimplemented\n";
    return SPoint2(0,0);
  }
}

int GFace::classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const
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
  }
  if(containsPoint(pt))
    return 1;
  return 0;
}

int GFace::classifyParam(const SPoint2 &par, int chkb, GEntity **bent) const
{
  if(chkb){
    AOMD::SPoint3 pt = point(par);
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
  }
  if(containsParam(par))
    return 1;
  return 0;
}

double * GFace::nthDerivative(const SPoint2 &, int,  double *) const
{ 
  cerr << "GFace::nthDerivative - not implemented\n"; 
  return nullptr;
}


GFPoint GFace::closestPoint(const AOMD::SPoint3 &)
{
  cerr << "GFace::closestPoint not implemented\n";
  return GFPoint();
}

int GFace::whichUse(GFaceUse const *use)
{ 
  if(use == &FrontUse)
    return 1;
  else if(use == &BackUse)
    return 0;
  cerr << "GFace::whichUse - error use not on face\n";
  return -1;
}

int GFace::inClosure(GEntity *ent) const
{
  return (this==ent || FrontUse.inClosure(ent) || BackUse.inClosure(ent));
}

SSList<GEdge*> GFace::splitU(double u)
{
  notImplemented();
  return SSList<GEdge*>();
}

SSList<GEdge*> GFace::splitV(double v)
{
  notImplemented();
  return SSList<GEdge*>();
}

SSList<GEdge*> GFace::splitSeamU(double ulow, double uhigh)
{
  notImplemented();
  return SSList<GEdge*>();
}

SSList<GEdge*> GFace::splitSeamV(double vlow, double vhigh)
{
  notImplemented();
  return SSList<GEdge*>();
}


void GFace::insertEdge(GEdge *e)
  // insert a new edge into the definition of a face, the new edge
  // points to the two vertices at its ends which must be on the 
  // closure of the face
  // this may split the face, returns a newly created face if this happens
{
  cerr << "Error GFace::insertEdge called\n";
}

void GFace::addLoop(SSList<GEdge*> edges, SSList<int> dirs)
{
  // need to make two loop uses
  GLoopUse *l1 = new GLoopUse(model(),&FrontUse,edges,dirs);
  FrontUse.addLoopUse(l1);
  edges.reverse();
  dirs.reverse();
  int d;
  for(SSListIter<int> dIter(dirs); dIter(d); )
    dIter.replaceCurrent(d?0:1);
  GLoopUse *l2 = new GLoopUse(model(),&BackUse,edges,dirs);
  BackUse.addLoopUse(l2);
  l1->setMate(l2);
  l2->setMate(l1);
}

void GFace::write(ostream &out)
{
  out << tag() << " " << id() << "\n  ";
  FrontUse.write(out);
  out << "\n  ";
  BackUse.write(out);
  out << "\n";
}

void GFace::read(istream &in)
{
  int id;
  in >> id;
  setID(id);
  FrontUse.read(in);
  BackUse.read(in);
}

SVector3 GFace::inwardSecant(GEdge *e, double ept, int dir)
  // this is just needed when building the topology of the model
  // to get the radial ordering of the faces around an edge

  // this is just returning the inward tangent to the face
  // really want a secant in the case where we have tangencies
{
  SPoint2 pt = e->reparamOnFace(this,ept,dir);
  SVector3 n = normal(pt);
  SVector3 etan = e->firstDer(ept);
  etan *= (dir ? 1.0 : -1.0);
  SVector3 ftan = cross(n,etan);
  ftan.normalize();
  return ftan;
}

double GFace::period(int dir) const
{
  if(!surfPeriodic(dir))
    return 0;
  else
    cerr << "GFace::period - error\n";
  return 0;
}

void GFace::fixPeriodicPar(SPoint2 &pt) const
  // fix the parameter in pt so that it is in the range of the face
  // if the underlying surface is periodic in that direction (shouldn't 
  // need to do anything otherwise. The reason for doing this if the
  // underlying surface is periodic, rather than the face, is that
  // the non-periodic face may lie on a periodic boundary of the surface
  // and some modelers will return that for a parameter point on the face
{
#ifdef DEBUG
  SPoint2 orig = pt;
#endif 

  // only trying to fix parameters that are off by a period
  // anything else means something else is wrong.
  
  for(int i=0; i < 2; i++){
    if(surfPeriodic(i)){
      // use this as a tolerance in parameter space
      AOMD::Range<double> range = parBounds(i);
      double tol = 1e-6*(range.high()-range.low());
      if(pt[i] < range.low()-tol)
	pt[i] += period(i);
      if(pt[i] > range.high()+tol)
	pt[i] -= period(i);
#ifdef DEBUG
      if(pt[i] < range.low()-tol){
	cerr << "GFace::fixPeriodicParam - error low " << i << "\n";
	cerr << "was: " << orig << " period = " << period(i);
	cerr << "fixed to: " << pt << "\n";
      }
      if(pt[i] > range.high()+tol){
	cerr << "GFace::fixPeriodicParam - error high " << i << "\n";
	cerr << "was: " << orig << " period = " << period(i);
	cerr << "fixed to: " << pt << "\n";
      }
#endif
      if(pt[i] < range.low())
	pt[i] = range.low();
      if(pt[i] > range.high())
	pt[i] = range.high();
    }
  }
}

void GFace::info() const
{
  cerr << name() << " - " << (void*)this << "\n";
  GFace *face = (GFace*)this;
  pGFaceUse fu = GF_use(face,1);
  pGFUIter fuIter = GFU_loopIter(fu);
  pGLoopUse lu;
  while(GFUIter_next(fuIter,&lu)){
    pGLUIter luIter = GLU_edgeIter(lu);
    pGEdgeUsePair e;
    int dir;
    cerr << "\n";
    while(GLUIter_next(luIter,&e,&dir))
      cerr << e->edge()->tag() << " " << dir << "\n";
  }
  cerr << "tag = " << face->tag() << "\n";
  double tol;
  C_modtol(2,face,&tol);
  cerr << "  tolerance = " << tol << "\n";
  int num;
  double par[2];
  C_parType2(2,face,0,&num,par);
  cerr << "orient = " << C_parOrient(Gface,face) << "\n";
  int ie = 0;
  GEdge *edge;
  for(SSListCIter<GEdge*> eIter(face->edges()); eIter(edge); ie++){
    cerr << "    edge " << ie << ": tag = " << edge->tag() << " ";
    cerr << "# faces: " << edge->faces().size() << "dir=";
    cerr << C_F_edgeDir(face,edge) << "\n";
    if(GE_isSeam(edge,face))
      cerr << "  edge is seam\n";
  }
  
  SSList<GVertex*> ve = face->vertices();
  //cerr << "  " << ve.size() << " vertices\n";
  GVertex *v;
  for(SSListCIter<GVertex*> vIter(ve); vIter(v); ){
    AOMD::SPoint3 point = v->point();
    //cerr << "  tag = " << v->tag() << "\n    ";
    //cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
  }
  SSList<GEdge*> vf = face->edges();
  cerr << "  " << vf.size() << " edges\n";
  SSList<GRegion*> vr = edge->regions();
  cerr << "  " << vr.size() << " regions\n";
  
  AOMD::Range<double> range = face->parBounds(0);
  cerr << "u range: " << range.low() << ", " << range.high() << "  ";
  cerr << "periodic = " << face->periodic(0) << "\n";
  double midx = 0.5*(range.low()+range.high());
  range = face->parBounds(1);
  cerr << "v range: " << range.low() << ", " << range.high() << "  ";
  cerr << "periodic = " << face->periodic(1) << "\n";
  double midy = 0.5*(range.low()+range.high());
  
  SPoint2 mid(midx,midy);
  SVector3 n = face->normal(mid);
  cerr << n << "\n";
  cerr << face->region(0) << " " << face->region(1) << "\n";
    
}




