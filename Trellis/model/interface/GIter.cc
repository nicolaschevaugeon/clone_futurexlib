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
#include "modeler.h"
//#include "GFace.h"
//#include "GVertex.h"
//#include "GFPoint.h"
//#include "GVPoint.h"
//#include "GEdge.h"
#include "SSList.h"
//#include "toPList.h"
//#include "Range.h"
//#include "SVector3.h"
////#include "SSList.h"
#include "SGModel.h"

extern "C"
GRIter GM_regionIter(pGModel mod)
{ return new SSListCIter<pGRegion>(mod->firstRegion());}

extern "C"
void GRIter_reset(GRIter iter)
{ iter->reset(); }

extern "C"
pGRegion GRIter_next(GRIter iter)
{
  pGRegion reg;
  int notDone = (*iter)(reg);
  if (notDone)
    return reg;
  return nullptr;
}

extern "C"
void GRIter_delete(GRIter iter)
{ delete iter; }

extern "C"
GFIter GM_faceIter(pGModel mod)
{ return new SSListCIter<pGFace>(mod->firstFace());}

extern "C"
void GFIter_reset(GFIter iter)
{ iter->reset(); }

extern "C"
pGFace GFIter_next(GFIter iter)
{
  pGFace face;
  int notDone = (*iter)(face);
  if (notDone)
    return face;
  return nullptr;
}

extern "C"
void GFIter_delete(GFIter iter)
{ delete iter; }

extern "C"
GEIter GM_edgeIter(pGModel mod)
{ return new SSListCIter<pGEdge>(mod->firstEdge());}

extern "C"
void GEIter_reset(GEIter iter)
{ iter->reset(); }

extern "C"
pGEdge GEIter_next(GEIter iter)
{
  pGEdge edge;
  int notDone = (*iter)(edge);
  if (notDone)
    return edge;
  return nullptr;
}

extern "C"
void GEIter_delete(GEIter iter)
{ delete iter; }

extern "C"
GVIter GM_vertexIter(pGModel mod)
{ return new SSListCIter<pGVertex>(mod->firstVertex());}

extern "C"
void GVIter_reset(GVIter iter)
{ iter->reset(); }

extern "C"
pGVertex GVIter_next(GVIter iter)
{
  pGVertex vtx;
  int notDone = (*iter)(vtx);
  if (notDone)
    return vtx;
  return nullptr;
}

extern "C"
void GVIter_delete(GVIter iter)
{ delete iter; }







