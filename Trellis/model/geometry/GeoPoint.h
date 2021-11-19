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
#ifndef H_GeoPoint
#define H_GeoPoint

/** A renderable representation of a model vertex as a point. */
class GeoPoint {
public:
  GeoPoint(double x, double y, double z);

  double x();
  double y();
  double z();

protected:
  double Location[3];
};

inline double GeoPoint::x()
{ return Location[0]; }

inline double GeoPoint::y()
{ return Location[1]; }

inline double GeoPoint::z()
{ return Location[2]; }

#endif
