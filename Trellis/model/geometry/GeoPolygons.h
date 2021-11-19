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
#ifndef H_GeoPolygons
#define H_GeoPolygons

/** A renderable representation of a model face represented as
   a bunch of polygons. */
class GeoPolygons {
public:
  GeoPolygons(int numPts);
  GeoPolygons(int numPts, double *pts,int nc, int *connect, int sizeP,
	      int numNorm, double *normals, int owns=0);
  ~GeoPolygons();

  int numPoints();
  int numPolys();
  int polySize(int i);
  int * poly(int i);
  double *point(int i);
  double *normal(int i); //returns normal at point i
protected:
  int Owns;
  int NumPts;
  int NumPolys;
  int PolySize;
  double *Pts;
  int *Connect;
  double *Normals;
};

inline int GeoPolygons::numPoints()
{ return NumPts; }

inline double * GeoPolygons::point(int i)
{ return &(Pts[i*3]); }

inline double * GeoPolygons::normal(int i)
{ return &(Normals[i*3]); }

inline int GeoPolygons::numPolys()
{ return NumPolys; }

inline int GeoPolygons::polySize(int)
{ return PolySize; }

inline int * GeoPolygons::poly(int i)
{ return &(Connect[i*PolySize]); }

#endif
