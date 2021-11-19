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
#ifndef H_GeoPolyLine
#define H_GeoPolyLine

/** A renderable representation of a model edge as a polyline. */
class GeoPolyLine {
public:
  GeoPolyLine(int numPts);
  GeoPolyLine(int numPts, double *pts);

  int numPoints();
  double *point(int i);
protected:
  int Owns;
  int NumPts;
  double *Pts;
};

inline int GeoPolyLine::numPoints()
{ return NumPts; }

inline double * GeoPolyLine::point(int i)
{ return &(Pts[i*3]); }

#endif
