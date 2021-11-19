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
// mVector.h: interface for the mVector class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MVECTOR_H_
#define _MVECTOR_H_

#include <cmath>
#include "mPoint.h"

#ifdef MSC
// #include "Constants.h"
#define M_PI 3.14159265358979323846
#endif

namespace Trellis_Util {
  class mTensor2;
  /**
     Simple class for a Vector
  */
  class mVector  
  {
    double pos[3];
  public:
    inline mVector(double x = 0.0, double y=0.0, double z=0.0)
    {pos[0]=x;pos[1]=y;pos[2]=z;}
    inline mVector(const mPoint &from, const mPoint &to);

    inline mVector operator- () const;        

    //jfadd const added
    inline mVector operator- (const mVector &other) const;
    //jfadd const added
    inline mVector operator+ (const mVector &other) const;
    inline mVector & operator += (const mVector &other);
    inline mVector & operator *= (const double &other);
    inline mVector & operator /= (const double &other);
    inline mVector & operator -= (const mVector &other);
    //jfadd const added
    inline double operator * (const mVector &other) const
    {
      return pos[0] * other.pos[0] + pos[1] * other.pos[1] + pos[2]*other.pos[2];
    }

    // LEFT mutiply of a mVector to a mTensor2. Added by T. Bui 2/03
    mVector operator * ( const mTensor2& M ) const;

    //jfadd const added
    inline mVector operator % (const mVector &other) const;
    //jfadd const added
    inline mVector operator * (const double &other) const;
    mVector & operator *= (const mTensor2 &other);
    //jfadd const added
    inline mVector operator / (const double &other) const;
    inline double operator () (int i) const;
    inline double & operator () (int i);
    inline double & operator [] (int i);
    inline double operator [] (int i) const;
    inline double angleRad (const mVector &v) const;
    inline double angleDeg (const mVector &v) const;
    //jfadd
    inline mVector& norm();
    inline double normValue();

    inline double mag() const;              // added by T. Bui 3/03

    inline ~mVector()= default;;
    friend class mTensor2;    
  };

  inline double mVector::mag() const
  {
          return sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
  }

  //jfadd
  inline double  mVector::normValue()
   {
    double n = ::sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    if(n==0.0)return n;
    pos[0]/=n;pos[1]/=n;pos[2]/=n;
    return n;
  }

  inline mVector&  mVector::norm()
   {
    double n = ::sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    if(n==0.0)return *this;
    pos[0]/=n;pos[1]/=n;pos[2]/=n;
    return *this;
  }

   inline double & mVector::operator () (int i)
  {
#ifdef _DEBUG_
    if(i>=3)throw new mException (__LINE__,__FILE__,"wrong index");
#endif
    return pos[i];
  }

  inline double mVector::operator () (int i) const
  {
    return pos[i];
  }

  inline double & mVector::operator [] (int i)
  {
#ifdef _DEBUG_
    if(i>=3)throw new mException (__LINE__,__FILE__,"wrong index");
#endif
    return pos[i];	    
  }

  inline double mVector::operator [] (int i) const
  {
#ifdef _DEBUG_
    if(i>=3)throw new mException (__LINE__,__FILE__,"wrong index");
#endif
    return pos[i];	    
  }

  inline mVector::mVector(const mPoint &from, const mPoint &to)
  {
    pos[0] = to(0) - from(0);
    pos[1] = to(1) - from(1);
    pos[2] = to(2) - from(2);
  }

  inline mVector mVector::operator - () const
  {
  //    added by T. Bui 1/03
                return mVector(-pos[0], -pos[1], -pos[2]);
  }

  //jfadd const added
  inline mVector mVector::operator - (const mVector &other) const
  {
    return mVector(pos[0]-other.pos[0],pos[1]-other.pos[1],pos[2]-other.pos[2]);
  }

  // cross product
  //jfadd const added
  inline mVector mVector::operator % (const mVector &other) const
  {
    return mVector(	pos[1]*other.pos[2]-pos[2]*other.pos[1],
			pos[2]*other.pos[0]-pos[0]*other.pos[2],
			pos[0]*other.pos[1]-pos[1]*other.pos[0]);
  }

  //jfadd const added
  inline mVector mVector::operator + (const mVector &other) const
  {
    return mVector(pos[0]+other.pos[0],pos[1]+other.pos[1],pos[2]+other.pos[2]);
  }

  //jfadd const added
  inline mVector mVector::operator * (const double &other) const
  {
    return mVector(pos[0]*other,pos[1]*other,pos[2]*other);
  }

  //jfadd const added
  inline mVector mVector::operator / (const double &other) const
  {
    return mVector(pos[0]/other,pos[1]/other,pos[2]/other);
  }

  inline mVector & mVector::operator += (const mVector &other)
  {
    pos[0]+=other.pos[0];pos[1]+=other.pos[1];pos[2]+=other.pos[2];
    return *this;
  }

  inline mVector & mVector::operator /= (const double &other)
  {
    double d = 1./other;
    pos[0]*=d;
    pos[1]*=d;
    pos[2]*=d;
    return *this;
  }

  inline mVector & mVector::operator -= (const mVector &other)
  {
    pos[0]-=other.pos[0];pos[1]-=other.pos[1];pos[2]-=other.pos[2];
    return *this;
  }

  inline mVector & mVector::operator *= (const double &other)
  {
    pos[0]*=other;pos[1]*=other;pos[2]*=other;
    return *this;
  }

  inline double mVector::angleRad(const mVector &v) const
  {
    // have to copy not to modify original vectors
    mVector y(v);
    mVector x(*this);
    x.norm();
    y.norm();
    double cosA = x * y;
    mVector cross = x % y;
    double sinA = ::sqrt(cross.pos[0]*cross.pos[0]+
			 cross.pos[1]*cross.pos[1]+
			 cross.pos[2]*cross.pos[2]);
    return atan2(sinA,cosA);
  }

  inline double mVector::angleDeg(const mVector &v) const
  {
     // return angleRad(v) * (180./3.14159265358979323846);
    return angleRad(v) * (180./M_PI);
  }

  
} // end of namespace

#endif 

