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
#ifndef H_SVector3
#define H_SVector3

#include "SPoint3.h"

/* concrete class for vector of size 3 */
class SVector3 {
public:
  SVector3() {}
  /// Construct from 2 SPoints, vector from p1 to p2
  SVector3(const AOMD::SPoint3 &p1, const AOMD::SPoint3 &p2)
    : P(p2-p1) {}
  /// Construct from a single SPoint, vector from origin to p1
  SVector3(const AOMD::SPoint3 &p1)
    : P(p1) {}
  SVector3(double x, double y, double z)
    : P(x,y,z) {}
  SVector3(double v)
    : P(v,v,v) {}
  SVector3(const double *array)
    : P(array) {}

  ///
  inline double x() const { return P.x(); }
  ///
  inline double y() const { return P.y(); }
  ///
  inline double z() const { return P.z(); }

  double normalize();

  // why both [] and (), why not
  double &operator[](int);
  double operator[](int) const;
  double &operator()(int);
  double operator()(int) const;

  ///
  SVector3 & operator += (const SVector3 &);
  ///
  SVector3 & operator -= (const SVector3 &);
  ///
  SVector3 & operator *= (const SVector3 &);
  ///
  SVector3 & operator = (double);

  operator double *() { return P; }
  
  friend std::istream& operator>>(std::istream& is, SVector3 &v);
  friend std::ostream& operator<<(std::ostream& os, SVector3 &v);

protected:
  AOMD::SPoint3 P;

};

///
double norm(const SVector3 &v);
///
double dot(const SVector3 &a, const SVector3 &b);
///
SVector3 cross(const SVector3 &a, const SVector3 &b);
///
double angle(const SVector3 &a, const SVector3 &b, const SVector3 &n);
///
SVector3 operator*(double m,const SVector3 &v);
///
SVector3 operator*(const SVector3 &v, double m);
///
SVector3 operator*(const SVector3 &v1, const SVector3 &v2);
///
SVector3 operator+(const SVector3 &a,const SVector3 &b);
///
SVector3 operator-(const SVector3 &a,const SVector3 &b);

inline double &SVector3::operator[](int i)
{ return P[i]; }

inline double SVector3::operator[](int i) const
{ return P[i]; }

inline double &SVector3::operator()(int i)
{ return P[i]; }

inline double SVector3::operator()(int i) const
{ return P[i]; }


#endif
