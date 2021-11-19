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
#include "TopoModel.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
#include "GeoRep.h"
#include "GeoPolygons.h"
#include "GVPoint.h"
#include "SPoint3.h"
#include "SSList.h"
#include "TEdge.h"
#include "TFace.h"

void main(int argc, char **argv)
{
  int i;
  
  if(argc != 2){
    cerr << "you forgot the filename!!!\n";
    exit(1);
  }
  SGModel *omodel = new ShapesModel(argv[1],0);
  omodel->setGeomTolerance(1.0);
  omodel->stats();
  //omodel->writeSMD();
  TopoModel *model = new TopoModel(omodel);

  GVertex *vtx;
  SPoint3 point;

  GFace *face;

  for(i=0 ; i < model->numFace(); i++){
    face = model->face(i);
    cerr << "face "<<i<<": tag = " << face->tag() << "\n";
    int periodicX = 0, periodicY=0;
    if(face->periodic(0) == Logical::True){
      cerr << "  periodic in u\n";
      periodicX = 1;
      Range<double> xrange = face->parBounds(0);
      cerr << "  u range from " << xrange.low() << " to " << xrange.high() << "\n";
      double middle = 0.5*(xrange.low()+xrange.high());
      face->splitU(middle);
    }
    else if(face->periodic(1) == Logical::True){
      cerr << "  periodic in v\n";
      periodicY = 1;
      Range<double> yrange = face->parBounds(1);
      cerr << "  v range from " << yrange.low() << " to " << yrange.high() << "\n";
      double middle = 0.5*(yrange.low()+yrange.high());
      face->splitV(middle);
    }
  }
  //face = model->face(0);
  //((TFace*)face)->splitU(0.0);
  
  GEdge *edge;
  //edge = model->edge(0);
  //double middle = 0.5*(edge->parBounds(0).low()+edge->parBounds(0).high());
  //edge->split(middle);

  exit(0);

  for(i=0 ; i < model->numVertex(); i++){
    vtx = model->vertex(i);
    
    point = vtx->point();
    cerr << "vertex "<<i<<": tag = " << vtx->tag() << "\n";
    cerr << "  number of uses = " << vtx->numUses() << "\n";
    cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
    
  }

  Range<double> range;

  for(i=0 ; i < model->numEdge(); i++){
    edge = model->edge(i);
    cerr << "edge "<<i<<": tag = " << edge->tag() << "\n";
    cerr << "  number of uses = " << edge->numUses() << "\n";
    for(int j = 0; j < 2; j++){
      vtx = edge->vertex(j);
      point = vtx->point();
      cerr << "    vertex "<<j<<": tag = " << vtx->tag() << "\n    ";
      cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
    }
    GeoRep *gr;
    gr = edge->geometry();
    if(gr){
      GeoPolyLine *gpl;
      gpl = gr->polyLine();
    }
  }


  for(i=0 ; i < model->numFace(); i++){
    face = model->face(i);
    cerr << "face "<<i<<": tag = " << face->tag() << "\n";
    int ie = 0;
    for(SSListCIter<GEdge*> eIter(face->edges()); eIter(edge); ie++)
    cerr << "    edge " << ie << ": tag = " << edge->tag() << "\n";

    Range<double> ub = face->parBounds(0);
    Range<double> vb = face->parBounds(1);
    cerr << "u range: " << ub.low() << " to " << ub.high() << "\n";
    cerr << "v range: " << vb.low() << " to " << vb.high() << "\n";

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
