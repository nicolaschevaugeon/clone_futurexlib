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

#ifndef H_GPoint
#define H_GPoint

#include "SPoint3.h"

namespace AOMD{
class GEntity;
}

/** Base class for points evaluated on geometric entities. */
class GPoint : public AOMD::SPoint3 {
public:

  virtual const AOMD::GEntity *gEnt() const = 0;
  int isValid() const;
  virtual GPoint *clone() const =0;

  // check if this id gets used anywhere.
  void setID(int id);
  int id();
protected:
  GPoint() {}
  GPoint(double x, double y, double z)
    : SPoint3(x,y,z) {}
  GPoint(const SPoint3 &pt)
    : SPoint3(pt) {}
  
private:
  int ID;
};

inline void GPoint::setID(int id)
{ ID = id; }

inline int GPoint::id()
{ return ID; }

inline int GPoint::isValid() const
{ return gEnt() != nullptr; }

#endif

