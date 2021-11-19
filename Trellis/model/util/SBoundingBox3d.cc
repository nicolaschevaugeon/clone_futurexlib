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
#include "SBoundingBox3d.h"
#include <cfloat>


namespace AOMD {

SBoundingBox3d::SBoundingBox3d()
: MinPt(DBL_MAX,DBL_MAX,DBL_MAX), MaxPt(-DBL_MAX,-DBL_MAX,-DBL_MAX)
{}

SBoundingBox3d::SBoundingBox3d(const AOMD::SPoint3 & pt)
: MinPt(pt), MaxPt(pt)
{}
void SBoundingBox3d::operator+=(const AOMD::SPoint3 &pt)
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

  if(pt[2] < MinPt[2])
    MinPt[2] = pt[2];
  if (pt[2] > MaxPt[2])
    MaxPt[2] = pt[2];

}

void SBoundingBox3d::operator+=(const SBoundingBox3d &box)
{
  (*this)+=box.MinPt;
  (*this)+=box.MaxPt;
}

AOMD::SPoint3 SBoundingBox3d::min() const
{ return MinPt; }

AOMD::SPoint3 SBoundingBox3d::max() const
{ return MaxPt; }

AOMD::SPoint3 SBoundingBox3d::center() const
{ return (MinPt+MaxPt)*.5; }

void SBoundingBox3d::operator*=(double scale)
{
  AOMD::SPoint3 center = (MinPt + MaxPt)*.5;
  MaxPt -= center;
  MinPt -= center;
  MaxPt *= scale;
  MinPt *= scale;
  MaxPt += center;
  MinPt += center;

}

void SBoundingBox3d::scale(double sx, double sy, double sz)
{
  AOMD::SPoint3 center = (MinPt + MaxPt)*.5;
  MaxPt -= center;
  MinPt -= center;
  MaxPt[0] *= sx;   MaxPt[1] *= sy; MaxPt[2] *= sz;
  MinPt[0] *= sx;   MinPt[1] *= sy; MinPt[2] *= sz;
  MaxPt += center;
  MinPt += center;
}

void SBoundingBox3d::makeCube()
{
  AOMD::SPoint3 len = MaxPt-MinPt;
  double scales[3];
  double max=-1.0;

  for(int i = 0; i < 3; i++)
    max = len[i] > max ? len[i] : max;

  for(int j = 0; j < 3; j++)
    scales[j] = max/len[j];

  scale(scales[0],scales[1],scales[2]);
  
}


}
