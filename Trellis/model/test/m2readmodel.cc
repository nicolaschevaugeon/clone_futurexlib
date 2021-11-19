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
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include "m2-interface/M2IModel.h"
#include "entities/GVertex.h"
#include "entities/GEdge.h"
#include "entities/GLoop.h"
#include "entities/GFace.h"
#include "geometry/GeoRep.h"
#include "geometry/GeoPolygons.h"
#include "entities/GVPoint.h"
#include "SPoint3.h"
#include "SSList.h"

void main(int argc, char **argv)
{
  int i;
  
  if(argc != 2){
    cout << "you forgot the filename!!!\n";
    exit(1);
  }
  SGModel *model = new M2IModel(argv[1]);
  model->stats();
  model->writeSMD();

  GVertex vtx;
  SPoint3 point;

  for(i=0 ; i < model->nVertex(); i++){
    vtx = model->vertex(i);
    
    point = vtx.point();
    cout << "vertex "<<i<<": tag = " << vtx.tag() << "\n";
    cout << point.x() << ", " << point.y() << ", " << point.z() << "\n";
    
  }

  Range<double> range;
  GEdge edge;

  for(i=0 ; i < model->nEdge(); i++){
    edge = model->edge(i);
    cout << "edge "<<i<<": tag = " << edge.tag() << "\n";
    //GeoRep *gr;
    //gr = edge.geometry();
    //GeoPolyLine *gpl;
    //gpl = gr->polyLine();
  }

  GFace face;

  for(i=0 ; i < model->nFace(); i++){
    face = model->face(i);
    cout << "face "<<i<<": tag = " << face.tag() << "\n";
    //GeoRep *gr;
    
    //gr = face.geometry();
    //GeoPolygons *gp;
    //gp = gr->polygon();
    //cerr << "Face " << i << "\n";
    //cerr << "Number of polygons: " << gp->numPolys() << "\n";

//     face = model->face(i);
//     SSList<GLoop> loops = face.loops();
//     cout << "Number of loops = " << loops.size() << "\n";

//     GLoop loop;
//     for(SSListIter<GLoop> lIter(loops); lIter(loop); )
//       cout << "  Loop: " << loop.numEdges() << " edges\n";

//     range = face.parBounds(0);
//     cout << range.low() << ", " << range.high() << "; ";

//     range = face.parBounds(1);
//     cout << range.low() << ", " << range.high() << "\n";

  }
return;
}
