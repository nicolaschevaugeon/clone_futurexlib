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
#include "SVector3.h"
#include <cmath>


SVector3 cross(const SVector3 &a, const SVector3 &b)
{
  return SVector3(a.y()*b.z()-b.y()*a.z(),
		 -(a.x()*b.z()-b.x()*a.z()),
		 a.x()*b.y()-b.x()*a.y());
}

double dot(const SVector3 &a, const SVector3 &b)
{
  return( a.x()*b.x()+a.y()*b.y()+a.z()*b.z() );
}

double angle(const SVector3 &a, const SVector3 &b, const SVector3 &n)
{
  double ddot = dot(a,b)/(norm(a)*norm(b));
  if(ddot>1.0)
    ddot = 1.0;
  if(ddot<-1.0)
    ddot = -1.0;
  double dang = acos(ddot);
  // check if we just found the angle or 2*PI-angle
  SVector3 check = cross(a,b);
  // ***** check this below ******
  if( norm(check) > 0.0000001*(norm(a)+norm(b))){  // check if a,b are colinear
    double dir = dot(check,n);
    dang = dir < 0 ? 2*M_PI-dang : dang;
  } else {
    if( ddot > 0)
      dang = 2*M_PI;
  }
  return dang;

}

double norm(const SVector3 &v)
{ return sqrt(dot(v,v)); }

double SVector3::normalize()
{ 
  double n = norm(*this);
  P[0] /= n; P[1]/= n; P[2] /= n;
  return n;
}

SVector3 operator*(double m,const SVector3 &v)
{
  return SVector3(v[0]*m,v[1]*m,v[2]*m);
}

SVector3 operator*(const SVector3 &v, double m)
{
  return SVector3(v[0]*m,v[1]*m,v[2]*m);
}

SVector3 operator*(const SVector3 &v1, const SVector3 &v2)
{
  return SVector3(v1[0]*v2[0],v1[1]*v2[1],v1[2]*v2[2]);
}

SVector3 operator+(const SVector3 &a,const SVector3 &b)
{
  return SVector3(a[0]+b[0],a[1]+b[1],a[2]+b[2]);
}

SVector3 operator-(const SVector3 &a,const SVector3 &b)
{
  return SVector3(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
}

SVector3 & SVector3::operator += (const SVector3 &a)
{ 
  for(int i = 0; i < 3; i++)
    P[i] += a[i];
  return *this;
}

SVector3 & SVector3::operator -= (const SVector3 &a)
{ 
  for(int i = 0; i < 3; i++)
    P[i] -= a[i];
  return *this;
}

SVector3 & SVector3::operator *= (const SVector3 &a)
{ 
  for(int i = 0; i < 3; i++)
    P[i] *= a[i];
  return *this;
}

SVector3 & SVector3::operator = (double v)
{
  P[0] = v; P[1] = v; P[2] = v;
  return *this;
}

istream& operator>>(istream& is, SVector3 &v)
{
  is >> v[0] >> v[1] >> v[2];
  return is;
}

ostream& operator<<(ostream& os, SVector3 &v)
{
  os << "[" << v[0] << ","<< v[1] << "," << v[2] << "]";
  return os;
}
