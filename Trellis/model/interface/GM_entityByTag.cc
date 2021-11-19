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
using std::cerr;
extern "C"
AOMD::GEntity * GM_entityByTag(SGModel *m,int type, int tag)
{
  switch(type){
  case TopoType::Vertex:
    return m->vertexByTag(tag);
  case TopoType::Edge:
    return m->edgeByTag(tag);
  case TopoType::Face:
    return m->faceByTag(tag);
  case TopoType::Region:
    return m->regionByTag(tag);
  }
  cerr << "GM_entityByTag - unknown entity type: "<< type << "\n";
  return nullptr;

}


