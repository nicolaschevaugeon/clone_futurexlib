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
#include "SBoundingBox2d.h"
#include <cfloat>

//#ifdef _SCOREC_NewCompiler
//namespace SCOREC_Util {
//#endif

SBoundingBox2d::SBoundingBox2d()
: MinPt(DBL_MAX,DBL_MAX), MaxPt(-DBL_MAX,-DBL_MAX)
{}

void SBoundingBox2d::operator+=(const SPoint2 &pt)
// note: it is possible for pt[i] to be both > MaxPt[i] and < MinPt[i]
// the first point always will be both
{
  if(pt[0] < MinPt[0])
    MinPt[0] = pt[0];
  if (pt[0] > MaxPt[0])
    MaxPt[0] = pt[0];

  if(pt[1] < MinPt[1])
    MinPt[1] = pt[1];
  if (pt[1] > MaxPt[1])
    MaxPt[1] = pt[1];

}

SPoint2 SBoundingBox2d::min() const
{ return MinPt; }

SPoint2 SBoundingBox2d::max() const
{ return MaxPt; }

SPoint2 SBoundingBox2d::center() const
{ return (MinPt+MaxPt)*.5; }

void SBoundingBox2d::operator*=(double scale)
{
  SPoint2 center = (MinPt + MaxPt)*.5;
  MaxPt -= center;
  MinPt -= center;
  MaxPt *= scale;
  MinPt *= scale;
  MaxPt += center;
  MinPt += center;

}

void SBoundingBox2d::scale(double sx, double sy)
{
  SPoint2 center = (MinPt + MaxPt)*.5;
  MaxPt -= center;
  MinPt -= center;
  MaxPt[0] *= sx;   MaxPt[1] *= sy;
  MinPt[0] *= sx;   MinPt[1] *= sy;
  MaxPt += center;
  MinPt += center;
}

void SBoundingBox2d::makeSquare()
{
  double xlen = MaxPt[0]-MinPt[0];
  double ylen = MaxPt[1]-MinPt[1];

  double ratio = xlen/ylen;
  if(ratio>1)
    scale(1.0,ratio);
  else
    scale(1/ratio,1.0);
  
}

//#ifdef _SCOREC_NewCompiler
//} // end of namespace SCOREC_Util
//#endif

