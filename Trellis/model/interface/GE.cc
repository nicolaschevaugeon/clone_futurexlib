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
#include "GFace.h"
#include "GVertex.h"
#include "GLoopUse.h"
#include "GEdgeUse.h"
#include "GEPoint.h"
#include "GVPoint.h"
#include "SSList.h"
#include "toPList.h"

extern "C"
int GE_periodic(GEdge *e)
{ return e->periodic(); }

extern "C"
int GE_inClosure(GEdge *e, AOMD::GEntity *ent)
{ return e->inClosure(ent); }

extern "C"
pPList GE_regions(GEdge *e)
{
  SSList<GRegion*> rs = e->regions();
  return toPList(rs);
}

extern "C"
pPList GE_faces(GEdge *e)
{
  SSList<GFace*> fs = e->faces();
  return toPList(fs);
}

extern "C"
pPList GE_vertices(GEdge *e)
{
  SSList<GVertex*> vs = e->vertices();
  return toPList(vs);
}

extern "C"
GVertex * GE_vertex(GEdge *e, int dir)
{
  return e->vertex(dir);
}

extern "C"
void GE_parRange(GEdge *e, double *low, double *high)
{
  AOMD::Range<double> b = e->parBounds(0);
  *low = b.low();
  *high = b.high();
}

extern "C"
int GE_point(GEdge *e, double par, double *xyz)
{
  GEPoint ep = e->point(par);
  ep.position(xyz);
  return ep.isValid();
}

extern "C"
double GE_vertexReparam(GEdge *e, GVertex *v)
{
  GVPoint vp(v->point(),v);
  return e->param(vp);
}

extern "C"
int GE_isSeam(GEdge *inedg, GFace *face)
{
  return inedg->isSeam(face);
}

extern "C"
pPList GE_uses(GEdge *e)
{ return toPList( e->uses() ); }

// these ops are unsupported and inefficient

extern "C"
GFace * GE_face(GEdge *e, int i)
{ return e->faces().nth(i); }

extern "C"
int GE_numFaces(GEdge *e)
{ return e->faces().size(); }


