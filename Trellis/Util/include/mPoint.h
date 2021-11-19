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
// mPoint.h: interface for the mPoint class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MPOINT_H_
#define _MPOINT_H_
#include <ostream>
namespace Trellis_Util {

/**
   Simple class for 3D points
 */

class mPoint  
{
  /// Position
  double pos[3];
 public:
  inline mPoint(double x = 0.0, double y =0.0, double z = 0.0)
    {
      pos[0] = x;
      pos[1] = y;
      pos[2] = z;
    }
  inline mPoint& operator +=(const mPoint &);
  inline mPoint operator +(const mPoint &) const;
  inline mPoint operator -(const mPoint &) const;
  inline mPoint operator *(double) const;
  inline mPoint& operator *=(double);
  inline ~mPoint() = default;
  inline double & operator () (int);
  inline double operator () (int) const;
  inline bool operator < (const mPoint &other) const;
  /** Lexicographic sorting i.e. sort by x then y then z.
      Cannot be used for geometrical search bu only for
      putting points in an associative container. Question
      like what is the closest point ot that point needs
  */
  inline bool lexicographicLessThan (const mPoint &other, double EPS) const;
};

inline std::ostream& operator<<(std::ostream& s, const Trellis_Util::mPoint& p)
{
  s << p(0) << " " << p(1) << " " << p(2);
  return s;
}

inline double & mPoint::operator () (int i) 
{
#ifdef _DEBUG_
	if(i>=3)throw new mException (__LINE__,__FILE__);
#endif
	return pos[i];
}

inline double mPoint::operator () (int i) const
{
#ifdef _DEBUG_
	if(i>=3)throw new mException (__LINE__,__FILE__);
#endif
	return pos[i];
}

inline mPoint mPoint::operator * (double other) const 
{
	return mPoint(pos[0]*other, pos[1]*other, pos[2]*other);
}

inline mPoint mPoint::operator + (const mPoint &other) const 
{
	return mPoint(pos[0]+other.pos[0], pos[1]+other.pos[1], pos[2]+other.pos[2]);
}

inline mPoint& mPoint::operator += (const mPoint &other) 
{
	pos[0]+=other.pos[0]; 
	pos[1]+=other.pos[1];
	pos[2]+=other.pos[2];
	return *this;
}

inline mPoint& mPoint::operator *= (double other) 
{
	pos[0]*=other; 
	pos[1]*=other;
	pos[2]*=other;
	return *this;
}

inline mPoint mPoint::operator - (const mPoint &other) const 
{
	return mPoint(pos[0]-other.pos[0], pos[1]-other.pos[1], pos[2]-other.pos[2]);
}

inline bool mPoint ::  lexicographicLessThan (const mPoint &other, double EPS) const
{
  if(pos[0]<other(0)-EPS)return 1;
  if(pos[0]>other(0)+EPS)return 0;
  if(pos[1]<other(1)-EPS)return 1;
  if(pos[1]>other(1)+EPS)return 0;
  if(pos[2]<other(2)-EPS)return 1;
  if(pos[2]>other(2)+EPS)return 0;
  return 0;
}

inline bool mPoint :: operator< (const mPoint &other) const
{
  return lexicographicLessThan(other,1.e-6);
}

} // end of namespace
#endif


