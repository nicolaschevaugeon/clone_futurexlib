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
#include "GVertex.h"
#include "GVPoint.h"
#include "SSList.h"
#include "toPList.h"

extern "C"
pPList GV_regions(GVertex *v)
{
  SSList<GRegion*> rs = v->regions();
  return toPList(rs);
}

extern "C"
pPList GV_faces(GVertex *v)
{
  SSList<GFace*> fs = v->faces();
  return toPList(fs);
}

extern "C"
pPList GV_edges(GVertex *v)
{
  SSList<GEdge*> es = v->edges();
  return toPList(es);
}

extern "C"
void GV_point(GVertex *v, double *xyz)
{
  GVPoint vp = v->point();
  vp.position(xyz);
}

extern "C"
pPList GV_uses(GVertex *e)
{ return toPList( e->uses() ); }

// these operators are inefficient and unsupported
extern "C"
GEdge * GV_edge(GVertex *v, int i)
{ return v->edges().nth(i); }

extern "C"
int GV_numEdges(GVertex *v)
{ return v->edges().size(); }


