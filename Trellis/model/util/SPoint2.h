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

#ifndef H_SPoint2
#define H_SPoint2

#include <cmath>
#include "iofwd.h"

/** A point in 2-space */
class SPoint2 {
protected:
  double P[2];
public:
  SPoint2() = default;
  ///
  SPoint2(double x, double y) 
    {P[0] = x; P[1] = y;}
  ///
  SPoint2(double *p)
    {P[0] = p[0]; P[1] = p[1];}
  SPoint2(const SPoint2 &pt)
    {P[0] = pt.P[0]; P[1] = pt.P[1]; }
  virtual ~SPoint2() = default;
  ///
  void setPosition(double xx, double yy);
  ///
  void getPosition(double *xx, double *yy) const;
  ///
  void position(double *) const;
  ///
  inline double x() const;
  ///
  inline double y() const;
  
  ///
  double &operator[](int);
  ///
  double operator[](int) const;
  SPoint2 &operator=(const SPoint2 &p);

  ///
  void operator+=(const SPoint2 &p);
  ///
  void operator-=(const SPoint2 &p);
  ///
  void operator*=(double mult);
  ///
  SPoint2 operator*(double mult);
  
  operator double *() { return P; }
};

///
double dist(const SPoint2 &, const SPoint2 &);
///
std::ostream & operator<<(std::ostream &out, const SPoint2 &p);

inline SPoint2 operator + (const SPoint2 &a, const SPoint2 &b)
{ return SPoint2(a.x()+b.x(),a.y()+b.y()); }

inline SPoint2 operator - (const SPoint2 &a, const SPoint2 &b)
{ return SPoint2(a.x()-b.x(),a.y()-b.y()); }

inline void SPoint2::setPosition(double xx, double yy)
{ P[0] = xx;  P[1] = yy; }

inline void SPoint2::getPosition(double *xx, double *yy) const
{ *xx = P[0];  *yy = P[1]; }

inline void SPoint2::position(double *p) const
{ p[0] = P[0]; p[1] = P[1]; }

inline double SPoint2::x() const
{ return P[0]; }

inline double SPoint2::y() const
{ return P[1]; }

inline SPoint2 & SPoint2::operator=(const SPoint2 &p)
{ P[0] = p.P[0]; P[1]=p.P[1]; return *this; }

inline double &SPoint2::operator[](int i)
{ return P[i]; }

inline double SPoint2::operator[](int i) const
{ return P[i]; }

inline void SPoint2::operator+=(const SPoint2 &p)
{ P[0] += p.P[0]; P[1] += p.P[1];}

inline void SPoint2::operator-=(const SPoint2 &p)
{ P[0] -= p.P[0]; P[1] -= p.P[1];}

inline void SPoint2::operator*=(double mult)
{ P[0] *= mult; P[1] *= mult; }

inline SPoint2 SPoint2::operator*(double mult)
{ return SPoint2(P[0]*mult, P[1]*mult); }


#endif

