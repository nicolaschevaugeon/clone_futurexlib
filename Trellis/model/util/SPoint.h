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

#ifndef H_SPoint
#define H_SPoint

#include <math.h>
#include <iostream>

template<int Dimension>
class SPoint {
protected:
  double P[Dimension];
public:
  SPoint() {}
  SPoint(double x)
    {P[0] = x;}
  SPoint(double x, double y)
    {P[0] = x; P[1] = y; }
  SPoint(double x, double y, double z) 
    {P[0] = x; P[1] = y; P[2] = z;}
  SPoint(double *p)
    {for(int i = 0; i < Dimension; i++) P[i] = p[i]; }
  SPoint(const SPoint<Dimension> &pt)
    {for(int i = 0; i < Dimension; i++) P[i] = pt.P[i]; }

  double &operator[](int);
  double operator[](int) const;
  SPoint<Dimension> &operator=(const SPoint<Dimension> &p);
    
  friend istream & operator>>(istream &in, SPoint<Dimension> &p);
};

template<int Dimension>
SPoint<Dimension>::SPoint(double *p)
{ for(int i = 0; i < Dimension; i++) P[i] = p[i]; }

template<int Dimension>
SPoint<Dimension>::SPoint(const SPoint<Dimension> &pt)
    {for(int i = 0; i < Dimension; i++) P[i] = pt.P[i]; }

template<int Dimension>
SPoint<Dimension> & SPoint<Dimension>::operator=(const SPoint<Dimension> &p)
{ 
  for(int i = 0; i < Dimension; i++) P[i] = p.P[i]; 
  return *this;
}

template<int Dimension>
double &SPoint<Dimension>::operator[](int i)
{ return P[i]; }

template<int Dimension>
double SPoint<Dimension>::operator[](int i) const
{ return P[i]; }


#endif

