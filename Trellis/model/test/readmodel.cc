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
#include "GEdgeUsePair.h"
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
  SGModel *model = new ShapesModel(argv[1],0);
  //model->setGeomTolerance(1.0);
  model->stats();
  model->writeSMD("test.smd");

  GVertex *vtx;
  GVPoint point;

  for(i=0 ; i < model->numVertex(); i++){
    vtx = model->vertex(i);
    
    point = vtx->point();
    //cerr << "vertex "<<i<<": tag = " << vtx->tag() << "\n";
    //cerr << "  number of uses = " << vtx->numUses() << "\n";
    //cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
    
    SSList<GEdge*> ve = vtx->edges();
    //cerr << "  " << ve.size() << " edges\n";
    SSList<GFace*> vf = vtx->faces();
    //cerr << "  " << vf.size() << " faces\n";
    SSList<GRegion*> vr = vtx->regions();
    //cerr << "  " << vr.size() << " regions\n";

  }
  Range<double> range;
  GEdge *edge;

  for(i=0 ; i < model->numEdge(); i++){
    edge = model->edge(i);
    cerr << "edge "<<i<<": tag = " << edge->tag() << "\n";
    //cerr << "  number of uses = " << edge->numUses() << "\n";
    for(int j = 0; j < 2; j++){
      vtx = edge->vertex(j);
      point = vtx->point();
      cerr << "    vertex "<<j<<": tag = " << vtx->tag() << "\n    ";
      cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
      cerr << "    param = " << edge->param(point) << "\n";
    }

    SSList<GVertex*> ve = edge->vertices();
    //cerr << "  " << ve.size() << " vertices\n";
    SSList<GFace*> vf = edge->faces();
    //cerr << "  " << vf.size() << " faces\n";
    SSList<GRegion*> vr = edge->regions();
    //cerr << "  " << vr.size() << " regions\n";
      
    GeoRep *gr;
    gr = edge->geometry();
    GeoPolyLine *gpl;
    gpl = gr->polyLine();
  }


  GFace *face;

  for(i=0 ; i < model->numFace(); i++){
    face = model->face(i);

    pGFaceUse fu = GF_use(face,1);
    pGFUIter fuIter = GFU_loopIter(fu);
    pGLoopUse lu;
    while(GFUIter_next(fuIter,&lu)){
      pGLUIter luIter = GLU_edgeIter(lu);
      pGEdgeUsePair e;
      int dir;
      cerr << "\n";
      while(GLUIter_next(luIter,&e,&dir))
	cerr << e->edge()->tag() << " " << dir << "\n";
    }
    cerr << "face "<<i<<": tag = " << face->tag() << "\n";
    double tol;
    C_modtol(2,face,&tol);
    cerr << "  tolerance = " << tol << "\n";
    int num;
    double par[2];
    C_parType2(2,face,0,&num,par);
    cerr << C_parOrient(Gface,face) << "\n";
    int ie = 0;
    for(SSListCIter<GEdge*> eIter(face->edges()); eIter(edge); ie++){
      cerr << "    edge " << ie << ": tag = " << edge->tag() << "\n";
      cerr << C_F_edgeDir(face,edge) << "\n";
      if(GE_isSeam(edge,face))
	cerr << "  edge is seam\n";
    }

    SSList<GVertex*> ve = face->vertices();
    //cerr << "  " << ve.size() << " vertices\n";
    GVertex *v;
    for(SSListCIter<GVertex*> vIter(ve); vIter(v); ){
      point = v->point();
      //cerr << "  tag = " << v->tag() << "\n    ";
      //cerr << point.x() << ", " << point.y() << ", " << point.z() << "\n";
    }
      
    SSList<GEdge*> vf = face->edges();
    cerr << "  " << vf.size() << " edges\n";
    SSList<GRegion*> vr = edge->regions();
    cerr << "  " << vr.size() << " regions\n";

#if 0    
    GeoRep *gr;    
    gr = face->geometry();
    GeoPolygons *gp;
    gp = gr->polygon();
    cerr << "Face " << i << "\n";
    cerr << "Number of polygons: " << gp->numPolys() << "\n";
#endif

    range = face->parBounds(0);
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
