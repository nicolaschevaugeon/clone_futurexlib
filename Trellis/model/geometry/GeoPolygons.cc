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
#include "GeoPolygons.h"

GeoPolygons::GeoPolygons(int numPts) :Owns(1), NumPts(numPts)
{ 
  Pts = new double[numPts]; 
  Connect = nullptr;
  Normals = nullptr;
}

GeoPolygons::GeoPolygons(int numPts, double *pts,int nc, int *connect, 
int sizeP, int numNorm, double *normals, int owns)
  :Owns(owns), NumPts(numPts), NumPolys(nc/sizeP), PolySize(sizeP), Pts(pts), Connect(connect), Normals(normals)
{}

GeoPolygons::~GeoPolygons()
{ 
  if(Owns){
    delete [] Pts;
    delete [] Connect;
    delete [] Normals;
  }
}
