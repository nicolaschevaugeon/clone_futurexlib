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
#ifndef H_SVector2
#define H_SVector2

#include "SPoint2.h"
class SVector3;

/** concrete class for vector of size 2 */
class SVector2 {
public:
  SVector2() {}
  ///
  SVector2(const SPoint2 &p1, const SPoint2 &p2)
    : P(p2-p1) {}
  ///
  SVector2(const SPoint2 &p1)
    : P(p1) {}
  ///
  SVector2(double x, double y)
    : P(x,y) {}
  SVector2(double *array)
    : P(array) {}

  ///
  inline double x() const { return P.x(); }
  ///
  inline double y() const { return P.y(); }

  ///
  double normalize();
  ///
  double &operator[](int);
  ///
  double operator[](int) const;

  void operator*=(double m);
protected:
  SPoint2 P;

};

///
double norm(const SVector2 &v);
///
double dot(const SVector2 &a, const SVector2 &b);
///
SVector3 cross(const SVector2 &a, const SVector2 &b);
///
double angle(const SVector2 &a, const SVector2 &b);
///
SVector2 operator*(double m,const SVector2 &v);
///
SVector2 operator*(const SVector2 &v, double m);
///
SVector2 operator+(const SVector2 &a,const SVector2 &b);
///
SVector2 operator-(const SVector2 &a,const SVector2 &b);

inline double &SVector2::operator[](int i)
{ return P[i]; }

inline double SVector2::operator[](int i) const
{ return P[i]; }


#endif
