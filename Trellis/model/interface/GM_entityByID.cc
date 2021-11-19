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
#include "GEntity.h"

using std::cerr;
extern "C"
AOMD::GEntity * GM_entityByID(SGModel *m,int type, int index)
{
  static SGModel *lastModel=nullptr;
  static int lastType=0;
  static int lastID=0;
  static AOMD::GEntity *lastEnt=nullptr;
  
  if(lastID == index && lastType == type && lastModel==m)
    return lastEnt;
  
  switch(type){
  case TopoType::Vertex:
    lastEnt =  m->vertexByID(index);
    break;
  case TopoType::Edge:
    lastEnt = m->edgeByID(index);
    break;
  case TopoType::Face:
    lastEnt = m->faceByID(index);
    break;
  case TopoType::Region:
    lastEnt = m->regionByID(index);
    break;
  default:
    cerr << "GM_entity - unknown entity type: "<< type << "\n";
    return nullptr;
  }
  lastModel = m;
  lastType = type;
  lastID = index;
  return lastEnt;
}


