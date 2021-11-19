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
#ifndef H_SBoundingBox2d
#define H_SBoundingBox2d

#include "SPoint2.h"

/** 2-d bounding box */
class SBoundingBox2d {
public:
  ///
  SBoundingBox2d();

  ///
  void operator+=(const SPoint2 &pt);

  ///
  void operator*=(double scale);

  ///
  void scale(double, double);
  ///
  SPoint2 min() const;
  ///
  SPoint2 max() const;
  ///
  SPoint2 center() const;
  ///
  void makeSquare();
private:
  SPoint2 MinPt,MaxPt;

};



#endif
