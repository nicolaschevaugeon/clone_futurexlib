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
#include "TopoModel.h"


extern "C"
SGModel *GM_createTopoModel(SGModel *model)
{
  return new TopoModel((SGModel*)model);
}

extern "C"
int GM_numVertices(SGModel *model)
{ return model->numVertex(); }

extern "C"
int GM_numEdges(SGModel *model)
{ return model->numEdge(); }

extern "C"
int GM_numFaces(SGModel *model)
{ return model->numFace(); }

extern "C"
int GM_numRegions(SGModel *model)
{ return model->numRegion(); }

extern "C"
GRegion * GM_outerRegion(SGModel *model)
{ return nullptr; }

extern "C"
void GM_writeSMD(SGModel *model, char *name)
{ model->writeSMD(name); }

extern "C"
double GM_tolerance(SGModel *model)
{ return model->tolerance(); }

extern "C"
void GM_bounds(SGModel *m, double *min, double *max)
{
  AOMD::SBoundingBox3d bounds = m->bounds();
  bounds.max().position(max);
  bounds.min().position(min);
}






