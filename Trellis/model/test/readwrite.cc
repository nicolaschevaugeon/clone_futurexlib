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
#include "ShapesModel.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
#include "GeoRep.h"
#include "GeoPolygons.h"
#include "GVPoint.h"
#include "SPoint3.h"
#include "SSList.h"

void main(int argc, char **argv)
{
  int i;
  
  if(argc != 2){
    cerr << "you forgot the filename!!!\n";
    exit(1);
  }
  ShapesModel::registerModeler();
  SGModel *model1 = new ShapesModel(argv[1],0);
  model1->setGeomTolerance(1.0);
  model1->stats();
  model1->writeSMD("test.smd");

  SGModel * model = SGModel::createModel("test.smd");
  model->writeSMD("test2.smd");

  GVertex *vtx;
  SPoint3 point;

  for(i=0 ; i < model->numVertex(); i++){
    vtx = model->vertex(i);
    
    point = vtx->point();
    cerr << "vertex "<<i<<": tag = " << vtx->tag() << "\n";
    cerr << "  number of uses = " << vtx->numUses() << "\n";
    cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
    
  }

  Range<double> range;
  GEdge *edge;

  for(i=0 ; i < model->numEdge(); i++){
    edge = model->edge(i);
    cerr << "edge "<<i<<": tag = " << edge->tag() << "\n";
    cerr << "  number of uses = " << edge->numUses() << "\n";
    for(int j = 0; j < 2; j++){
      vtx = edge->vertex(j);
      point = vtx->point();
      cerr << "    vertex "<<j<<": tag = " << vtx->tag() << "\n    ";
      //cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
    }
      
    GeoRep *gr;
    gr = edge->geometry();
    GeoPolyLine *gpl;
    gpl = gr->polyLine();
  }


  GFace *face;

  for(i=0 ; i < model->numFace(); i++){
    face = model->face(i);
    cerr << "face "<<i<<": tag = " << face->tag() << "\n";
    int ie = 0;
    for(SSListCIter<GEdge*> eIter(face->edges()); eIter(edge); ie++)
    cerr << "    edge " << ie << ": tag = " << edge->tag() << "\n";
    
    GeoRep *gr;    
    gr = face->geometry();
    GeoPolygons *gp;
    gp = gr->polygon();
    cerr << "Face " << i << "\n";
    cerr << "Number of polygons: " << gp->numPolys() << "\n";

//     face = model->face(i);
//     SSList<GLoop> loops = face->loops();
//     cerr << "Number of loops = " << loops.size() << "\n";

//     GLoop loop;
//     for(SSListIter<GLoop> lIter(loops); lIter(loop); )
//       cerr << "  Loop: " << loop.numEdges() << " edges\n";

//     range = face->parBounds(0);
//     cerr << range.low() << ", " << range.high() << "; ";

//     range = face->parBounds(1);
//     cerr << range.low() << ", " << range.high() << "\n";

  }
return;
}
