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
#include "GVertex.h"
#include "GFPoint.h"
#include "GEPoint.h"
#include "GVPoint.h"
#include "GEdgeUse.h"
#include "GLoopUse.h"
#include "GEdge.h"
#include "SSList.h"
#include "toPList.h"


extern "C"
int GF_geomDirection(GFace *f)
{ return f->geomDirection();  }

extern "C"
int GF_paramDegeneracies(GFace *f, int dir, double *par)
{ return f->paramDegeneracies(dir,par); }

extern "C"
int GF_inClosure(GFace *f, AOMD::GEntity *ent)
{ return f->inClosure(ent); }

extern "C"
pPList GF_regions(GFace *f)
{
  SSList<GRegion*> rs = f->regions();
  return toPList(rs);
}

extern "C"
pPList GF_edges(GFace *f)
{
  SSList<GEdge*> es = f->edges();
  return toPList(es);
}

extern "C"
pPList GF_vertices(GFace *f)
{
  SSList<GVertex*> vs = f->vertices();
  return toPList(vs);
}

extern "C"
GFaceUse * GF_use(GFace *f, int dir)
{ return f->use(dir); }

extern "C"
GRegion * GF_region(GFace *f, int dir)
{ return f->region(dir); }

extern "C"
void GF_parRange(GFace *e, int axis, double *low, double *high)
{
  AOMD::Range<double> b = e->parBounds(axis);
  *low = b.low();
  *high = b.high();
}

extern "C"
int GF_point(GFace *f, double *par, double *xyz)
{
  GFPoint fp = f->point(par[0],par[1]);
  fp.position(xyz);
  return fp.isValid();
}

extern "C"
void GF_normal(GFace *f, double *par, double *xyz)
{
  SVector3 n = f->normal(SPoint2(par[0],par[1]));
  xyz[0] = n[0];
  xyz[1] = n[1];
  xyz[2] = n[2];
}

extern "C"
void GF_edgeReparam(GFace *f, GEdge *e, double epar, int edir, double *fpar)
{
  //GEPoint ep(e->point(epar),e,epar);
  SPoint2 fpt = e->reparamOnFace(f,epar,edir);
  fpar[0] = fpt[0];
  fpar[1] = fpt[1];
}

extern "C"
void GF_vertexReparam(GFace *f, GVertex *v, double *fpar)
{
  GVPoint vp(v->point(),v); // dumb to have to give point and vertex
  SPoint2 fpt = f->param(vp);
  fpar[0] = fpt[0];
  fpar[1] = fpt[1];

}

extern "C"
pPList GF_splitU(GFace *f, double u)
{
  SSList<GEdge*> es = f->splitU(u);
  return toPList(es);
}

extern "C"
pPList GF_splitV(GFace *f, double v)
{
  SSList<GEdge*> es = f->splitV(v);
  return toPList(es);
}

extern "C"
int GF_periodic(GFace *f, int dir)
{ return f->periodic(dir); }



