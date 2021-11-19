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
#include "GEdge.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GLoopUse.h"
#include "GFace.h"
#include "GShell.h"
#include "GRegion.h"
#include "GEdgeUsePair.h"
#include "SSList.h"
#include "SGModel.h"
#include "GEPoint.h"
#include "GVPoint.h"
#include "SVector3.h"
#include "MessageOut.h"

using std::cerr;

GEdge::GEdge(SGModel *model, int t, GVertex *v0, GVertex *v1)
: GEntity(model,t)
{
  Vertices[0] = v0;
  Vertices[1] = v1;

}

GEdge::~GEdge()
= default;

TopoType::Value GEdge::type() const
{ return TopoType::Edge; }

int GEdge::dim() const
{ return 1; }

SSList<GRegion*> GEdge::regions() const
{
  SSList<GRegion*> rs;
  GEdgeUsePair * eu;
  for(SSListCIter<GEdgeUsePair*> euIter(Uses);  euIter(eu); ){
    GRegion *r;
    if((r=eu->shell()->region()))
      rs.appendUnique(r);
  }
  return rs;
}

SSList<GFace*> GEdge::faces() const
{
  SSList<GFace*> fs;
  GEdgeUsePair * eu;
  for(SSListCIter<GEdgeUsePair*> euIter(Uses);  euIter(eu); ){
    for(int i =0; i < 2; i++){
      GFace *f = eu->loopUse(i)->faceUse()->face();
      fs.appendUnique(f);
    }
  }
  return fs;
}

SSList<GEdge*> GEdge::edges() const
{
  SSList<GEdge*> es;
  return es;
}

SSList<GVertex*> GEdge::vertices() const
{
  SSList<GVertex*> vs;
  vs.append(vertex(0));
  vs.append(vertex(1));
  return vs;
}

GVertex* GEdge::vertex( int n) const
{ return Vertices[n]; }

Logical::Value GEdge::periodic(int dim) const
{ 
  if(vertex(0) == vertex(1))
    return Logical::True;
  return Logical::False;
}

Logical::Value GEdge::continuous(int dim) const
{
  return Logical::True;
}

int GEdge::isSeam(GFace *face) const
  // this is really brute force and should be overridden in derived classes
  // if there is a better way, also relies on parFromEdgePar to behave
  // nicely if edge isn't being used in direction passed
{
  AOMD::Range<double> r = parBounds(0);
  double mid = 0.5*(r.low()+r.high());
  SPoint2 fpt0 = reparamOnFace(face,mid,0);
  SPoint2 fpt1 = reparamOnFace(face,mid,1);
  // note: the comparison is being done in parametric space, thus this
  // isn't the right tolerance to use, but it's the only one we've got
  if(dist(fpt0,fpt1) < tolerance())
    return 0;
  else 
    return 1;
}

double GEdge::param(const GPoint & pt)
{
  switch( pt.gEnt()->dim() ){
  case 0: {
    const GEntity *gent = pt.gEnt();
    AOMD::Range<double> edgeRange = parBounds(0);
    if(gent == vertex(0))
      return (geomDirection() == 1 ? edgeRange.low() : edgeRange.high());
    else if(gent == vertex(1))
      return (geomDirection() == 1 ? edgeRange.high() : edgeRange.low());
    else
      cerr << "Edge::param - vertex not on edge\n";
    return 0;
  }
  case 1:
    return ((GEPoint &)pt).par();
  default:
    cerr << "GEdge::param unimplemented for point on face or region\n";
    return 0;
  }
}

int GEdge::classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const
{
  GVertex *v;
  if(chkb){
    for (SSListCIter<GVertex *> vIter(vertices()); vIter(v); ){
      if(v->containsPoint(pt)){
	if(bent)
	  *bent = v;
	return -1;
      }
    }
  }
  if(containsPoint(pt))
    return 1;
  return 0;
}

int GEdge::classifyParam(double par, int chkb, GEntity **bent) const
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
  }
  if(containsParam(par))
    return 1;
  return 0;
}


GEPoint GEdge::closestPoint(const AOMD::SPoint3 &)
{
  cerr << "Edge::closestPoint not implemented\n";
  return GEPoint();
}

void GEdge::setVertex(GVertex *v, int which)
{ Vertices[which] = v; }

GEdgeUse * GEdge::addLoop(GLoopUse *l, int dir)
{ 
  GEdgeUsePair *eu;
  if(Uses.size() == 1 && !Uses.nth(0)->loopUse(0) &&
     !Uses.nth(0)->loopUse(1)){ // the first use was created with the edge, 
                                //no loop set, (I don't think this applies
                                // anymore
    eu = Uses.nth(0);
    eu->setLoopUse(l,dir);
  } else {
    eu = new GEdgeUsePair(model(),this,l,dir,Vertices[0],Vertices[1]);
    Uses.append(eu); 
  }
  return eu->side(dir);
}

int GEdge::inClosure(GEntity *ent) const
{
  return (this==ent || Vertices[0]==ent || Vertices[1]==ent);
}

void GEdge::removeUse(GEdgeUsePair *eu)
{
  Uses.remove(eu);
}

#if 0
void GEdge::mergeUses(GShell *su)
  // merge the uses of this edge w.r.t the given shell use
{
  const double angleTol = 0.0000001;  // not sure 

  GEdgeUsePair *eui;
  SSList<GEdgeUsePair*> toBeMerged;
  for(SSListCIter<GEdgeUsePair*> euIter(Uses); euIter(eui); ){
    if(eui->shell()==su)
      toBeMerged.append(eui);
  }
  if(toBeMerged.size() == 2)
    toBeMerged.nth(0)->mergeWith(toBeMerged.nth(1));
  else if(toBeMerged.size() == 4){
    // need to merge the right pairs, find dirs and face uses
    GEdgeUsePair *eu0;
    eu0 = toBeMerged.nth(0);
    int dir = eu0->loopUse(0) ? 0 : 1;
    if(!eu0->loopUse(dir))
      cerr << "GEdge::mergeUses - error 1\n";
    GFace *f = eu0->loopUse(dir)->faceUse()->face();
    int mw0;  // index of use to merge with use 0
    for(mw0=1; mw0 < 4; mw0++){
      eui = toBeMerged.nth(mw0);
      if( eui->loopUse(!dir) && eui->loopUse(!dir)->faceUse()->face() != f)
	break;
    }
    if(mw0 == 4){ 
      // this case is when a cylindrical face is used on both sides
      // by the same shell, merge with same loop use in other direction
      GLoopUse *lu0 = eu0->loopUse(dir);
      for(mw0=1; mw0 < 4; mw0++){
	eui = toBeMerged.nth(mw0);
	if( eui->loopUse(!dir) == lu0 )
	  break;
      }
      if(mw0==4)
	cerr << "GEdge::mergeUses - error 2\n";
    } 
    eu0->mergeWith(toBeMerged.nth(mw0));
    // now merge the two remaining
    GEdgeUsePair *eu2 = toBeMerged.nth( mw0 == 1 ? 2 : 1);
    GEdgeUsePair *eu3 = toBeMerged.nth( mw0 == 1 ? 3 : (mw0 == 2 ? 3 : 2) );
    // check some stuff, debugging
    int d2 = eu2->loopUse(0) ? 0 : 1;
    int d3 = eu3->loopUse(0) ? 0 : 1;
    if(d2 == d3)
      cerr << "GEdge::mergeUses - error 3\n";
    eu2->mergeWith(eu3);
  } else{  // not 2 or 4 uses
    int nuses = toBeMerged.size();
    GEdgeUsePair **eus = new GEdgeUsePair*[nuses+2];
    double *angles = new double[nuses/2+1];
    for(int i=0; i < nuses+2; i++)
      eus[i]=0; 
    for(int j=0; j < nuses/2+1; j++)
      angles[j]=0; 
    
    Range<double> erange = parBounds(0);
    double middle = 0.5*(erange.low()+erange.high());
    SVector3 etan = firstDer(middle);
    SVector3 ref;
    int first = 1;

    GEdgeUsePair *eu;
    int nu=0;
    for(SSListIter<GEdgeUsePair*> euIter(toBeMerged); euIter(eu); ){
      int dir = eu->loopUse(0) ? 0 : 1;  // direction use is used in
      if(!eu->loopUse(dir))
	cerr << "GEdge::mergeUses - error 11\n";
      GFaceUse *fu = eu->loopUse(dir)->faceUse();
      int fudir = fu->dir();
      int fdir = !(dir ^ fudir);  // dir face is using edge through this use
      GFace *f = fu->face();

      SVector3 fnormal = f->inwardSecant(this,middle,fdir);

      int same = 0;
      if(first){
	ref = fnormal;
	eus[dir] = eu;
	nu++;
	first = 0;
      } else {
	double a = angle(ref,fnormal,etan);
	if( fabs(a-2*M_PI) < angleTol)
	  a = 0.0;
	// insert into array sorted
	int iu;
	for(iu = 0; iu < nu; iu++){
	  if( fabs(angles[iu] - a) < angleTol ){
	    same = 1;
	    break;
	  }
	  if(angles[iu] > a)
	    break;
	}
	// this is a member of the iu'th pair
	if(!same){
	  nu++;
	  // first make room, remember array is nuses+2 in size
	  int imove;
	  for(imove = nuses-1; imove >= iu*2; imove--) // shift up 2
	    eus[imove+2] = eus[imove];
	  eus[iu*2] = 0;
	  eus[iu*2+1] = 0;
	  for(imove = nuses/2-1; imove >= iu; imove--) // shift up 1
	    angles[imove+1] = angles[imove];
	  angles[iu] = a;
	}
	if(eus[iu*2+dir])
	  cerr << "already have this use\n";
	eus[iu*2+dir] = eu;
      }
    }
    // have uses in order of increasing angle in the array
    // use i is paired with use (i+1)%(nuses+2), if use i == 0 then
    // its mate should be also
    
    // it's possible that the array got filled such that the last two
    // positions are zero, check this and adjust lastUse so that
    // the last and first use will still get paired up
    int lastUse = nuses+2;
    if( eus[nuses]==0 && eus[nuses+1]==0)
      lastUse = nuses;

    int iuse = 1;
    int usesRemaining = nuses/2;
    while(usesRemaining){
      int imate = (iuse+1)%lastUse;
      GEdgeUsePair *eu2=0, *eu3=0;
      eu2 = eus[iuse];
      eu3 = eus[imate];
      if(eu2==0 && eu3==0){
	iuse+=2;
	continue;
      }
      if(eu2 ==0 || eu3 ==0){
	cerr << "warning GEdge::mergeUses - single use\n";
	// something bad has happened here
	// assume that there is a match for this use and find it
	int ifind,izero;
	if(eu2==0){
	  izero = ifind = iuse;
	} else {
	  izero = ifind = imate;
	}
	while(ifind<lastUse && !eus[izero]){
	  ifind+=2;
	  if(eus[ifind]){  
	    eus[izero] = eus[ifind];
	    eus[ifind] = 0;
	    if(eu2==0)
	      eu2 = eus[izero];
	    else
	      eu3 = eus[izero];
	  }
	}
	if( eus[nuses]==0 && eus[nuses+1]==0)
	  lastUse = nuses;
      }
      // check some stuff, debugging
      int d2 = eu2->loopUse(0) ? 0 : 1;
      int d3 = eu3->loopUse(0) ? 0 : 1;
      if(d2 == d3)
	cerr << "GEdge::mergeUses - direction - error\n";
      eu2->mergeWith(eu3);
      iuse+=2;
      usesRemaining--;
    }
  }
}
#endif

void GEdge::mergeUses(GShell *su)
  // merge the uses of this edge 
{
  //const double angleTol = 0.0000001;  // not sure 

  GEdgeUsePair *eui;
  SSList<GEdgeUsePair*> toBeMerged;
  for(SSListCIter<GEdgeUsePair*> euIter(Uses); euIter(eui); ){
    if(eui->shell()==su && ((eui->loopUse(0)) ==nullptr || (eui->loopUse(1)==nullptr)))
      toBeMerged.append(eui);
  }
  if(toBeMerged.size() == 0)
    return;
  if(toBeMerged.size() == 2)
    toBeMerged.nth(0)->mergeWith(toBeMerged.nth(1));
  else if(toBeMerged.size() == 4){
    // need to merge the right pairs, find dirs and face uses
    GEdgeUsePair *eu0;
    eu0 = toBeMerged.nth(0);
    int dir = eu0->loopUse(0) ? 0 : 1;
    if(!eu0->loopUse(dir))
      cerr << "GEdge::mergeUses - error 1\n";
    GFace *f = eu0->loopUse(dir)->faceUse()->face();
    int mw0;  // index of use to merge with use 0
    for(mw0=1; mw0 < 4; mw0++){
      eui = toBeMerged.nth(mw0);
      if( eui->loopUse(!dir) && eui->loopUse(!dir)->faceUse()->face() != f)
	break;
    }
    if(mw0 == 4){ 
      // this case is when a cylindrical face is used on both sides
      // by the same shell, merge with same loop use in other direction
      GLoopUse *lu0 = eu0->loopUse(dir);
      for(mw0=1; mw0 < 4; mw0++){
	eui = toBeMerged.nth(mw0);
	if( eui->loopUse(!dir) == lu0 )
	  break;
      }
      if(mw0==4)
	cerr << "GEdge::mergeUses - error 2\n";
    } 
    eu0->mergeWith(toBeMerged.nth(mw0));
    // now merge the two remaining
    GEdgeUsePair *eu2 = toBeMerged.nth( mw0 == 1 ? 2 : 1);
    GEdgeUsePair *eu3 = toBeMerged.nth( mw0 == 1 ? 3 : (mw0 == 2 ? 3 : 2) );
    // check some stuff, debugging
    int d2 = eu2->loopUse(0) ? 0 : 1;
    int d3 = eu3->loopUse(0) ? 0 : 1;
    if(d2 == d3)
      cerr << "GEdge::mergeUses - error 3\n";
    eu2->mergeWith(eu3);
  } else{  // not 2 or 4 uses
    //cerr << "GEdge::mergeUses - in new part of code\n";
    //int nuses = toBeMerged.size();

    //int nfaces = nuses/2;
    SBlock<GFace*> faces;
    SBlock<int> dirs;

    int nf = getOrderedFaces(&faces,&dirs);
    // this returns an array of faces and a matching array of
    // directions with the face-direction pairs ordered correctly
    // around the edge

    while(toBeMerged.size()){
      GEdgeUsePair *eu = toBeMerged.nth(0);
      //int nu=0;
      int dir = eu->loopUse(0) ? 0 : 1;  // direction use is used in
      if(!eu->loopUse(dir))
	cerr << "GEdge::mergeUses - error 11\n";
      GFaceUse *fu = eu->loopUse(dir)->faceUse();
      int fudir = fu->dir();
      int fdir = !(dir ^ fudir);  // dir face is using edge through this use
      GFace *f = fu->face();
      
      // have face-direction pair, need to find in list
      int found = -1;
      for(int i = 0; i < nf; i++){
	if(faces[i]==f && dirs[i]==fdir){
	  found = i;
	  break;
	}
      }
      if(found == -1)
	cerr << "GEdge::mergeUses() - didn't find face\n";
      
      // based on which faceuse we're on we need to go either in
      // the positive direction or negative direction around the edge
      int rdir = (dir == 1 ? 1: nf-1);
      int mergeId = (found+rdir)%nf;
      GFace * mergeFace = faces[mergeId];
      int mergeDir = dirs[mergeId];

      // and the use we're looking for depends on the edge use directions
      int mergeUseDir;
      if(mergeDir == fdir)
	mergeUseDir = !fudir;
      else
	mergeUseDir = fudir;

      // need to find the use with mergeFace and mergeDir;
      GEdgeUsePair *meu, *mergeUse = nullptr;
      for(SSListIter<GEdgeUsePair*> meuIter(toBeMerged); meuIter(meu); ){
	int mdir = meu->loopUse(0) ? 0 : 1;  // direction use is used in
	if(!meu->loopUse(mdir))
	  cerr << "GEdge::mergeUses - error 11\n";
	GFaceUse *mfu = meu->loopUse(mdir)->faceUse();
	int mfudir = mfu->dir();
	int mfdir = !(mdir ^ mfudir);  // dir face is using edge with this use
	GFace *mf = mfu->face();
	if(mf==mergeFace && mfdir == mergeDir && mergeUseDir == mfudir){
	  mergeUse = meu;
	  break;
	}
      }
      if(!mergeUse)
	cerr << "GEdge::mergeUses - error: couldn't find merge use\n";
      
      eu->mergeWith(mergeUse);
      toBeMerged.remove(eu);
      toBeMerged.remove(mergeUse);
    }
  }
}

int GEdge::getOrderedFaces(SBlock<GFace*> *fs, SBlock<int> *dirs) const
{
  // can't be more than Uses.size faces at edge
  fs->reserve(Uses.size());
  dirs->reserve(Uses.size());

  SBlock<double> angles(Uses.size());

  AOMD::Range<double> erange = parBounds(0);
  double middle = 0.5*(erange.low()+erange.high());
  SVector3 etan = firstDer(middle);
  SVector3 ref;
  int nf = 0;

  GEdgeUsePair * eu;
  for(SSListCIter<GEdgeUsePair*> euIter(Uses); euIter(eu); ){
    int dir = eu->loopUse(0) ? 0 : 1;  // direction use is used in
    if(!eu->loopUse(dir))
      cerr << "GEdge::getOrderedFaces - error 11\n";
    GFaceUse *fu = eu->loopUse(dir)->faceUse();
    int fudir = fu->dir();
    int fdir = !(dir ^ fudir);  // dir face is using edge through this use
    GFace *f = fu->face();
    
    // check if we already have info for this face in this direction
    int have = 0;
    for(int i = 0; i < nf; i++){
      if((*fs)[i]==f && (*dirs)[i]==fdir){
	have = 1;
	break;
      }
    }
    if(have)
      continue;

    SVector3 fnormal = f->inwardSecant((GEdge*)this,middle,fdir);
    
    // int same = 0;
    if(nf==0){
      ref = fnormal;
      (*fs)[0] = f;
      (*dirs)[0] = fdir;
      angles[0] = 0;
    } else {
      double a = angle(ref,fnormal,etan);
      (*fs)[nf] = f;
      (*dirs)[nf] = fdir;
      angles[nf] = a;
    }
    nf++;
  }

  // sort fs,dirs on angles
  for(int i = 0; i < nf; i++){
    for(int j = 0; j < nf-i-1; j++){
      if(angles[j+1] < angles[j]){
	GFace *f = (*fs)[j+1];
	double a = angles[j+1];
	int dir = (*dirs)[j+1];
	(*fs)[j+1] = (*fs)[j];
	angles[j+1] = angles[j];
	(*dirs)[j+1] = (*dirs)[j];
	(*fs)[j] = f;
	angles[j] = a;
	(*dirs)[j] = dir;
      }
    }
  }
  return nf;
}

int GEdge::numUses() const
{ return Uses.size(); }

SSList<GEdgeUsePair*> GEdge::uses() const
{ return Uses; }

void GEdge::replace(GEdge *enew)
  // replace this edge with enew, updating all adjacency info
{
  enew->Uses = Uses;
  GEdgeUsePair *eu;
  for(SSListIter<GEdgeUsePair*> euIter(Uses); euIter(eu); )
    eu->replaceOwner(enew);
  // don't need to update vertices since they were given when edge
  // was created
}

void GEdge::splitAtVertex(GVertex *v, GEdge *e1, GEdge *e2)
{
  GEdgeUsePair *eu;
  for(SSListIter<GEdgeUsePair*> euIter(Uses); euIter(eu); )
    e2->Uses.append(eu->split(v,e2));
  replace(e1);  // replace this with e1
  model()->remove(this);
  delete this;

}

void GEdge::write(ostream &out)
{
  int i = 0;
  GEdgeUsePair *eu;
  out << tag() << " " << id() << " " << Uses.size() << " ";
  out << Vertices[0]->id() << " " << Vertices[1]->id() << "\n";
  for(SSListCIter<GEdgeUsePair*> euIter(Uses); euIter(eu); ){
    eu->setID(i++);
    eu->write(out);
  }
  out << "\n";
}

void GEdge::read(istream &in)
{
  int nuse, id;
  in>> id >> nuse;
  setID(id);
  
  int v0, v1;
  in >> v0 >> v1;
  Vertices[0] = model()->vertexByID(v0);
  Vertices[1] = model()->vertexByID(v1);

  for(int i = 0; i < nuse; i++){
    GEdgeUsePair *eu = new GEdgeUsePair(model(),this);
    Uses.append(eu);
    eu->setID(i);
    eu->read(in);
  }
}

GEdgeUsePair * GEdge::use(int id)
{ return Uses.nth(id); }

void GEdge::info() const
{
  cerr << name() << " - " << (void*)this << "\n";
  cerr << "vertices:\n";
  cerr << "  " << vertex(0)->name() << " - " << vertex(0) << "\n";
  cerr << "  " << vertex(1)->name() << " - " << vertex(1) << "\n";
  cerr << "uses:\n";
  GEdgeUsePair *eu;
  for(SSListCIter<GEdgeUsePair *> euIter(Uses); euIter(eu); )
    eu->info();
  
}

double GEdge::period() const
{
  if(!periodic(0))
    return 0;
  else
    cerr << "GEdge::period - error\n";
  return 0;
}

