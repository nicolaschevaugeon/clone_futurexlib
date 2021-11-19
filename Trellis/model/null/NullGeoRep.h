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

#ifndef H_NullGeoRep
#define H_NullGeoRep
#include "GeoRep.h"

class NullEntity;
class GeoPolygons;
class GeoPolyLine;
class GeoPoint;

/** This really can't do much of anything. */
class NullGeoRep : public GeoRep {
public:
  NullGeoRep( AOMD::GEntity *e, int d) : GeoRep(e), Dim(d) {}

  int topoDim() override;

  GeoPolygons * polygon() override;

  GeoPolyLine * polyLine() override;
  
  GeoPoint * point() override;

protected:
  int Dim;


};

#endif
