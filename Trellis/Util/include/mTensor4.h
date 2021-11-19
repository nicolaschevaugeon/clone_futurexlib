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
// mTensor4.h: interface for the mTensor4 class.
//
//////////////////////////////////////////////////////////////////////
#ifndef _MTensor4_H_
#define _MTensor4_H_

#include "mTensor2.h"

namespace Trellis_Util {

/**
   Simple class for a 4th order tensor in 3d which
   has a 3x3x3x3 matrix representation.
*/
  
  class mTensor4  
  {
    double pos[3][3][3][3];
  public:
    mTensor4 ()= default;
    mTensor4 (const double c[3][3][3][3]);
    mTensor4 (const double E, const double nu);
    inline double operator () (int i,int j, int k, int l) const;
    inline double & operator () (int i,int j, int k, int l);
    inline mTensor2 operator * (const mTensor2 &other) const
    {
      mTensor2 m(0.0);
      for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	  m(i,j) = 0.0;
	  for(int k=0;k<3;k++){
	    for(int l=0;l<3;l++){
	      m(i,j) += other(k,l) * pos[k][l][i][j];
	    }
	  }
	}
      }
      return m;
    }
  };
  
  inline double delta(int i, int j) {return ((i == j)?1.0:0.0);}

  inline mTensor4::mTensor4 (const double E, const double nu)
  {
    const double lam = nu*E/((1.+nu)*(1.-2.*nu));
    const double mu  = E/(2.*(1.+nu));
    int i, j, k, l;
    for (i = 0; i < 3; ++i){
      for (j = 0; j < 3; ++j){
	for (k = 0; k < 3; ++k){
	  for (l = 0; l < 3; ++l){
	    pos[i][j][k][l] = lam * delta(i,j) * delta(k,l) 
	      + mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	  }
	}
      }
    }
    return;
  }
  
  inline mTensor4::mTensor4 (const double c[3][3][3][3])
  {
    int i, j, k, l;
    for (i = 0; i < 3; ++i){
      for (j = 0; j < 3; ++j){
        for (k = 0; k < 3; ++k){
          for (l = 0; l < 3; ++l){
            pos[i][j][k][l] = c[i][j][k][l];
          }
        }
      }
    }
    return;
  }

  inline double mTensor4::operator () (int i,int j, int k, int l) const
  {
    return pos[i][j][k][l];
  }
  inline double & mTensor4::operator () (int i,int j, int k, int l)
  {
    return pos[i][j][k][l];
  }
} // end of namespace

#endif
