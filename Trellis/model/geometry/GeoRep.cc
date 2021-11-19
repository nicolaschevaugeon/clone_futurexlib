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
#include "GeoRep.h"
#include "GeoPolyLine.h"
#include "GeoPoint.h"
#include "SBlock.h"
#include "Range.h"
#include "SPoint3.h"
#include "GEdge.h"
#include "GVertex.h"
#include "GEPoint.h"
#include "GVPoint.h"
#include "GEdgeUtil.h"

GeoRep::~GeoRep()
= default;

int GeoRep::topoDim()
{ return Ent->dim(); }

GeoPolygons * GeoRep::polygon()
{ return nullptr; }

GeoPolyLine * GeoRep::polyLine()
{
  if(Ent->dim() != 1)
    return nullptr;

  GEdge *edge = (GEdge*)(Ent);
  //double shortest = 100000;
  AOMD::Range<double> range = edge->parBounds(0);

  SBlock<double> points(1200);
  int curPoint = 0;
  
  double ps,pe;
  ps = range.low();  pe = range.high();

  //  int count = 0;
  int done = 0;
  double d,pnl;
  
  AOMD::SPoint3 firstPoint = edge->point( range.low() );
  points[0] = firstPoint[0]; points[1] = firstPoint[1]; 
  points[2] = firstPoint[2];
  curPoint++;
 
  double middle = .5*(range.low()+range.high()); // middle of edge
  // calc the length of the edge (very roughly)
  double len = 2*dist(edge->point(range.low()),edge->point(middle));
  if(len > 0.0){ 
    
    d = len/20.0; 
    
    pnl = parAtDist(edge,pe,d,.5*d,-1);
    
    double pn, pl = ps;
    while(!done){
      int ok = 1;
      pn = parAtDist(edge,pl,d,.5*d,1);
      
      if(pn >= pnl)
	break;
      
      if(ok){
        AOMD::SPoint3 p = edge->point(pn);
	if(points.size() < (curPoint+1)*3)
	  points.reserve(curPoint*3*2);
	points[3*curPoint] = p[0]; 
	points[3*curPoint+1] = p[1]; 
	points[3*curPoint+2] = p[2];
	curPoint++;
	pl = pn;
      }
    }
  } else {
    AOMD::SPoint3 p = edge->vertex(0)->point();
    points[0] = p[0];
    points[1] = p[1];
    points[2] = p[2];
    curPoint = 1;
  }
  AOMD::SPoint3 p = edge->point( range.high() );
  points[3*curPoint] = p[0]; 
  points[3*curPoint+1] = p[1]; 
  points[3*curPoint+2] = p[2];
  curPoint++;

  double * pointList = new double[curPoint*3];
  for(int i = 0; i < curPoint*3; i++)
    pointList[i] = points[i];
  return new GeoPolyLine(curPoint,pointList);

}

GeoPoint * GeoRep::point()
{
  if(Ent->dim() != 0)
    return nullptr;
  
  GVertex *vertex = (GVertex*)Ent;
  AOMD::SPoint3 pt = vertex->point();
  return new GeoPoint(pt[0],pt[1],pt[2]);
}
