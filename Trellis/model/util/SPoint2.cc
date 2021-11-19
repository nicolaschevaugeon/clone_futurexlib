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
#include "SPoint2.h"

double dist(const SPoint2 &a, const SPoint2 &b)
{ return sqrt( (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])); }

ostream & operator<<(ostream &out, const SPoint2 &p)
{ return out << "(" << p.x() << "," << p.y() << ")"; }
