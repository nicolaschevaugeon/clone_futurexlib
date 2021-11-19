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
#include "SVector2.h"
#include "SVector3.h"
#include <cmath>

void SVector2::operator*=(double m)
{ P[0] *=m; P[1]*=m; }

SVector3 cross(const SVector2 &a, const SVector2 &b)
{
  return SVector3(0,0,a.x()*b.y()-b.x()*a.y());
}

double dot(const SVector2 &a, const SVector2 &b)
{
  return( a.x()*b.x()+a.y()*b.y() );
}

double angle(const SVector2 &a, const SVector2 &b)
{
  SVector3 n(0,0,1);
  double ddot = dot(a,b)/(norm(a)*norm(b));
  if(ddot>1.0)
    ddot = 1.0;
  if(ddot<-1.0)
    ddot = -1.0;
  double dang = acos(ddot);
  // check if we just found the angle or 2*PI-angle
  SVector3 check = cross(a,b);
  if( norm(check) > 0.01*(norm(a)+norm(b))){  // check if a,b are colinear
    double dir = dot(check,n);
    dang = dir < 0 ? 2*M_PI-dang : dang;
  } else {
    if( ddot > 0)
      dang = 2*M_PI;
  }
  return dang;

}

double norm(const SVector2 &v)
{ return sqrt(dot(v,v)); }

double SVector2::normalize()
{ 
  double n = norm(*this);
  P[0] /= n; P[1]/= n;
  return n;
}

SVector2 operator*(double m,const SVector2 &v)
{
  return SVector2(v[0]*m,v[1]*m);
}

SVector2 operator*(const SVector2 &v, double m)
{
  return SVector2(v[0]*m,v[1]*m);
}

SVector2 operator+(const SVector2 &a,const SVector2 &b)
{
  return SVector2(a[0]+b[0],a[1]+b[1]);
}

SVector2 operator-(const SVector2 &a,const SVector2 &b)
{
  return SVector2(a[0]-b[0],a[1]-b[1]);
}


