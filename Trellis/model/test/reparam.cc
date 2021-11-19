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
#include "modeler.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
#include "GFPoint.h"
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
  SGModel *model = new ShapesModel(argv[1],0);
  model->setGeomTolerance(1.0);
  model->stats();
  model->writeSMD("test.smd");

  GVertex *vtx;
  GVPoint point;

  GFace *face;

  for(i=0 ; i < model->numFace(); i++){
    face = model->face(i);
    cerr << "face "<<i<<": tag = " << face->tag() << "\n";
    SPoint3 chkpt = face->point(-3,2);
    cerr << " param(-3,2) is " << chkpt << "\n";
    chkpt = face->point(3,2);
    cerr << " param(3,2) is " << chkpt << "\n";

    cerr << C_parOrient(Gface,face) << "\n";
    int ie = 0;

    SSList<GEdge*> vf = face->edges();
    cerr << "  " << vf.size() << " edges\n";

    Range<double> range = face->parBounds(0);
    cerr << "u range: " << range.low() << ", " << range.high() << "  ";
    cerr << "periodic = " << face->periodic(0) << "\n";
    double midx = 0.5*(range.low()+range.high());
    range = face->parBounds(1);
    cerr << "v range: " << range.low() << ", " << range.high() << "  ";
    cerr << "periodic = " << face->periodic(1) << "\n";
    double midy = 0.5*(range.low()+range.high());

    SPoint2 mid(midx,midy);
    SVector3 n = face->normal(mid);
    cerr << n << "\n";
    cerr << face->region(0) << " " << face->region(1) << "\n";

    GEdge *edge;
    for(SSListCIter<GEdge*> eIter(face->edges()); eIter(edge); ie++){
      cerr << "    edge " << ie << ": tag = " << edge->tag() << "\n";
      cerr << C_F_edgeDir(face,edge) << "\n";
      Range<double> erange = edge->parBounds(0);
      double ecenter = 0.25*(erange.low()+erange.high());
      cerr << "    edge quarter parameter = " << ecenter << "\n";
      SPoint2 fpt;
      fpt = face->parFromEdgePar(edge,ecenter,0);
      cerr << "    dir = 0 : param = " << fpt << "\n";
      fpt = face->parFromEdgePar(edge,ecenter,1);
      cerr << "    dir = 1 : param = " << fpt << "\n";
    }


  }

  GRegion *region;
  for(i=0 ; i < model->numRegion(); i++){
    region = model->region(i);
    //cerr << "region "<<i<<":  = " << region << "\n";

    SSList<GVertex*> ve = region->vertices();
    //cerr << "  " << ve.size() << " vertices\n";
    SSList<GEdge*> vf = region->edges();
    //cerr << "  " << vf.size() << " edges\n";
    SSList<GFace*> vr = region->faces();
    //cerr << "  " << vr.size() << " faces\n";

  }
  
  return;
}
