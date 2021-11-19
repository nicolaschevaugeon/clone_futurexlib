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
// mTensor2.h: interface for the mTensor2 class.
//
//////////////////////////////////////////////////////////////////////
#ifndef _MTensor2_H_
#define _MTensor2_H_

#include "mVector.h"

namespace Trellis_Util
{
class mTensor4;
/**
   Simple class for a 2nd order tensor in 3d which
   has a 3x3 matrix representation.
*/

class mTensor2
{
   // double pos[3][3];
   mVector pos[3];

  public:
   /// Constuct a mTensor2, each line defined by the value inside vector
   mTensor2(const mVector &l1, const mVector &l2, const mVector &l3);
   /// construct a mTensor2 with all entries = init.
   mTensor2(double init = 0.0);
   mTensor2 operator*(const mTensor4 &other) const;

   inline mVector operator*(const mVector &other) const
   {
      mVector m(0, 0, 0);
      /*
       for(int i=0;i<3;i++)
         {
           for(int j=0;j<3;j++)
             {
               m.pos[i] += pos[i][j] * other.pos[j];
             }
         }
       */
      for (int i = 0; i < 3; ++i)
      {
         m[i] += pos[i] * other;
      }

      return m;
   }
   inline mTensor2 operator*(const mTensor2 &other) const
   {
      mTensor2 m(0);
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            for (int l = 0; l < 3; l++) m.pos[i][j] += pos[i][l] * other.pos[l][j];
         }
      }
      return m;
   }
   inline mTensor2 &operator*=(const double &scalar)
   {
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            pos[i][j] *= scalar;
         }
      }
      return *this;
   }

   inline mTensor2 &operator/=(const double &scalar)
   {
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            pos[i][j] /= scalar;
         }
      }
      return *this;
   }

   /*    inline mTensor2 operator = (const double &scalar) const
   {
     mTensor2 other(scalar);
     return other;
     }*/

   inline mTensor2 &operator+=(const mTensor2 &other)
   {
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            pos[i][j] += other.pos[i][j];
         }
      }
      return *this;
   }
   // jfadd
   inline mTensor2 &operator-=(const mTensor2 &other)
   {
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            pos[i][j] -= other.pos[i][j];
         }
      }
      return *this;
   }

   inline mTensor2 operator*(const double &scalar) const
   {
      mTensor2 other;
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            other.pos[i][j] = scalar * pos[i][j];
         }
      }
      return other;
   }
   double &operator()(int, int);
   double operator()(int, int) const;
   /// return line i of the tensor.
   const mVector &operator[](int i) const;
   mVector &operator[](int i);
   /// return col j of the tensor.
   mVector col(int j) const;

   mTensor2 operator+(const mTensor2 &other) const;
   virtual ~mTensor2();
   /// solve the system
   mVector operator/(const mVector &other) const;
   mTensor2 invert() const;
   /// determinant
   double det() const;
   double trace() const;
   double trace2() const;
   long eigen(mVector e[3], double v[3]) const;
   long eigen2d(mVector e[3], double v[3]) const;
   inline void symmetrize();
   mTensor2 intersect2d(mTensor2 &t) const;
   // mTensor2 transpose() const;

   void transpose();
   mTensor2 operator!() const;
};

inline void mTensor2::symmetrize()
{
   pos[0][1] = pos[1][0] = 0.5 * (pos[0][1] + pos[1][0]);
   pos[0][2] = pos[2][0] = 0.5 * (pos[0][2] + pos[2][0]);
   pos[2][1] = pos[1][2] = 0.5 * (pos[2][1] + pos[1][2]);
}

inline mTensor2::mTensor2(const mVector &c0, const mVector &c1, const mVector &c2)
{
   /*
       pos[0][0] = c0(0);
       pos[0][1] = c0(1);
       pos[0][2] = c0(2);
       pos[1][0] = c1(0);
       pos[1][1] = c1(1);
       pos[1][2] = c1(2);
       pos[2][0] = c2(0);
       pos[2][1] = c2(1);
       pos[2][2] = c2(2);
   */
   pos[0] = c0;
   pos[1] = c1;
   pos[2] = c2;
}

inline double &mTensor2::operator()(int i, int j) { return pos[i][j]; }

inline double mTensor2::operator()(int i, int j) const { return pos[i][j]; }

inline const mVector &mTensor2::operator[](int i) const
{
   //    return the ith row
   return pos[i];
}

inline mVector &mTensor2::operator[](int i)
{
   //    assign the ith row
   return pos[i];
}

inline mVector mTensor2::col(int j) const { return mVector(pos[0][j], pos[1][j], pos[2][j]); }

// NICO NEW FUNCTIONS
inline mTensor2 mTensor2::operator+(const mTensor2 &other) const
{
   mTensor2 m(0);
   /*
     for(int i=0;i<3;i++)
       {
         for(int j=0;j<3;j++)
           {
             m.pos[i][j] = pos[i][j] + other.pos[i][j];
           }
       }
    */
   for (int i = 0; i < 3; ++i) m.pos[i] += pos[i] + other[i];

   return m;
}
// NICO END
long FindCubicRoots(const double coeff[4], double x[3]);
long NullSpace(const double *a, double *result, double eps, long n);

inline mVector operator*(const mVector &other, const mTensor2 &t)
{
   mVector m(0, 0, 0);
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         m[i] += t(i, j) * other[j];
      }
   }
   return m;
}

}  // namespace Trellis_Util

#endif
