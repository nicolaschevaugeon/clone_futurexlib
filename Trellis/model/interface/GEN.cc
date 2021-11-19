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
#include "GEntity.h"
using std::cerr;

class GVertex;
class GEdge;
class GFace;

extern "C" {
void GV_point(GVertex* v, double *xyz);
double GE_vertexReparam(GEdge *e, GVertex *v);
int GE_point(GEdge* e, double par, double *xyz);
int GF_point(GFace* f, double *par, double *xyz);
void GF_edgeReparam(GFace *f, GEdge *e, double epar, int edir, double *fpar);
void GF_vertexReparam(GFace *f, GVertex *v, double *fpar);
}

extern "C"
double GEN_tolerance(AOMD::GEntity *e)
{
  return e->tolerance();
}

extern "C"
void GEN_bounds(AOMD::GEntity *e, double *min, double *max)
{
  AOMD::SBoundingBox3d bounds = e->bounds();
  bounds.max().position(max);
  bounds.min().position(min);
}

extern "C"
int GEN_repType(AOMD::GEntity *e)
{
  return e->repType();
}

extern "C"
void * GEN_getNativePtr(AOMD::GEntity *e)
{ return e->getNativePtr(); }

extern "C"
int GEN_getNativeInt(AOMD::GEntity *e)
{ return e->getNativeInt(); }

extern "C"
int GEN_inClosure(AOMD::GEntity *target, AOMD::GEntity *ent)
{ return target->inClosure(ent); }

extern "C"
void GEN_reparam(AOMD::GEntity *target, AOMD::GEntity *ent, double *entPar, int entDir,
		 double *targetPar) {

#ifdef DEBUG
  // Check to see if ent is in the closure of the target
  if (!GEN_inClosure(target, ent)) {
    cerr << "GEN_reparam - entity is not in the closure of the target entity\n";
    return;
  }
#endif

  // Get the type of target
  switch (target->type()) {
  case TopoType::Region:
    cerr << "GEN_reparam - can not reparam onto a region\n";
    return;
  case TopoType::Face:
    switch (ent->type()) {
    case TopoType::Edge :
      GF_edgeReparam((GFace *) target, (GEdge *) ent, *entPar, entDir, targetPar);
      return;
    case TopoType::Vertex :
      GF_vertexReparam((GFace *) target, (GVertex *) ent, targetPar);
      return;
    default:
      cerr << "GEN_reparam - face reparam failure\n";
      return;
    }
  case  TopoType::Edge:
    *targetPar = GE_vertexReparam((GEdge *) target, (GVertex *)ent);
      return;
  default:
    cerr << "GEN_reparam - reparam failure\n";
    return;
  }
  
}

extern "C"
int GEN_point(AOMD::GEntity *ent, double *entPar, double *xyz) {
  
  // Get the type of ent
  switch (ent->type()) {
  case TopoType::Region:
    cerr << "GEN_point - can not evaluate a point in a region\n";
    return 0;
  case TopoType::Face:
    return GF_point((GFace *) ent, entPar, xyz);
  case TopoType::Edge:
    return GE_point((GEdge *) ent, entPar[0], xyz);
  case TopoType::Vertex:
    GV_point((GVertex *) ent, xyz);
    return 1;
  default:
    cerr << "GEN_point - unsupported entity type\n";
    return 0;
  }
}

