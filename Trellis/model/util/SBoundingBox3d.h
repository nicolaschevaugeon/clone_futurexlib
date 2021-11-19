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
#ifndef H_SBoundingBox3d
#define H_SBoundingBox3d

#include "SPoint3.h"

namespace AOMD {
/** A bounding box class - add points and it grows to be the
  bounding box of the point set
  */
class SBoundingBox3d {
public:
  ///
  SBoundingBox3d();
  ///
  SBoundingBox3d(const AOMD::SPoint3 &);

  ///
  void operator+=(const AOMD::SPoint3 &pt);
  ///
  void operator+=(const SBoundingBox3d &pt);
  ///
  void operator*=(double scale);
  ///
  void scale(double, double, double);
  ///
  AOMD::SPoint3 min() const;
  ///
  AOMD::SPoint3 max() const;
  ///
  AOMD::SPoint3 center() const;
  ///
  void makeCube();
private:
  AOMD::SPoint3 MinPt,MaxPt;

};

}
#endif
