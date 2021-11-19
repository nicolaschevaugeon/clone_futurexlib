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
#ifndef H_GeoRep
#define H_GeoRep

class GeoPolygons;
class GeoPolyLine;
class GeoPoint;
namespace AOMD{
  class GEntity;
}

/** Encapsulates the displable representation of a model entity. 
 This is actually poorly named. */
class GeoRep {
public:
  GeoRep(AOMD::GEntity *ent) : Ent(ent) {}
  virtual ~GeoRep();

  /** Get the dimension of the entity. */
  virtual int topoDim();

  /** Get representation as polygons, if it exists. */
  virtual GeoPolygons * polygon();

  /** Get representation as polylines, if it exists. */
  virtual GeoPolyLine * polyLine();
  
  /** Get the representation as a point, if it exists. */
  virtual GeoPoint * point();

  /** Get the entity. */
  AOMD::GEntity *entity();
protected:
  AOMD::GEntity *Ent;


};

inline AOMD::GEntity * GeoRep::entity()
{ return Ent; }

#endif
