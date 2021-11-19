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
// mTensor2.cpp: implementation of the mTensor2 class.
//
//////////////////////////////////////////////////////////////////////

#include "mTensor2.h"
#include "mTensor4.h"
#include "mVector.h"
#include <algorithm>
#include <functional>
#include <cmath>
#include <cstdio>

#ifndef M_PI
#define M_PI 3.1415926
#endif

namespace Trellis_Util {

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

  mTensor2::mTensor2(double init)
  {
    for(int i=0;i<3;i++)
      {
	for(int j=0;j<3;j++)
	  {
	    pos[i][j] = init;
	  }
      }
  }
  
  mTensor2::~mTensor2()
  = default;
  
  double mTensor2::trace () const
  {
    return pos[0][0] + pos[1][1] + pos[2][2];
  }

  double mTensor2::trace2 () const
  {
    double a00 =  pos[0][0] * pos[0][0] + 
      pos[1][0] * pos[0][1] + 
      pos[2][0] * pos[0][2]; 
    double a11 =  pos[1][0] * pos[0][1] + 
      pos[1][1] * pos[1][1] + 
      pos[1][2] * pos[2][1]; 
    double a22 =  pos[2][0] * pos[0][2] + 
      pos[2][1] * pos[1][2] + 
      pos[2][2] * pos[2][2];

    return a00 + a11 + a22;
  }

  double mTensor2::det () const
  {
    return pos[0][0] * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]) -
      pos[0][1] * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]) +
      pos[0][2] * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
  }
  
  mVector mTensor2::operator / (const mVector & b) const
  {
    mVector res;
    
    double detm = det();
    
    if (detm == 0.0){
      throw;
    }
    
    double ud = 1. / (detm);
    
    res[0] = b[0] * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]) -
      pos[0][1] * (b[1] * pos[2][2] - pos[1][2] * b[2]) +
      pos[0][2] * (b[1] * pos[2][1] - pos[1][1] * b[2]);
    
    res[1] = pos[0][0] * (b[1] * pos[2][2] - pos[1][2] * b[2]) -
      b[0] * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]) +
      pos[0][2] * (pos[1][0] * b[2] - b[1] * pos[2][0]);
    
    res[2] = pos[0][0] * (pos[1][1] * b[2] - b[1] * pos[2][1]) -
      pos[0][1] * (pos[1][0] * b[2] - b[1] * pos[2][0]) +
      b[0] * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
    
    for (int i = 0; i < 3; i++)
      res[i] *= ud;
    
    return res;
  }

  mTensor2 mTensor2::invert() const
  {
    mTensor2 inv;
    
    double detm = det();
    
    if (detm == 0.0){
      throw;
    }
    
    double ud = 1. / (detm);
    inv.pos[0][0] = ud * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]);
    inv.pos[0][1] = -ud * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]);
    inv.pos[0][2] = ud * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
    inv.pos[1][0] = -ud * (pos[0][1] * pos[2][2] - pos[0][2] * pos[2][1]);
    inv.pos[1][1] = ud * (pos[0][0] * pos[2][2] - pos[0][2] * pos[2][0]);
    inv.pos[1][2] = -ud * (pos[0][0] * pos[2][1] - pos[0][1] * pos[2][0]);
    inv.pos[2][0] = ud * (pos[0][1] * pos[1][2] - pos[0][2] * pos[1][1]);
    inv.pos[2][1] = -ud * (pos[0][0] * pos[1][2] - pos[0][2] * pos[1][0]);
    inv.pos[2][2] = ud * (pos[0][0] * pos[1][1] - pos[0][1] * pos[1][0]);
    return inv;
  }


  
  
  long FindCubicRoots(const double coeff[4], double x[3])
  {
    double a1 = coeff[2] / coeff[3];
    double a2 = coeff[1] / coeff[3];
    double a3 = coeff[0] / coeff[3];
    
    double Q = (a1 * a1 - 3 * a2) / 9.;
    double R = (2. * a1 * a1 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
    double Qcubed = Q * Q * Q;
    double d = Qcubed - R * R;

    //    printf ("d = %22.15e Q = %12.5E R = %12.5E Qcubed %12.5E\n",d,Q,R,Qcubed);

    /// three roots, 2 equal 
    if(Qcubed == 0.0 || fabs ( Qcubed - R * R ) < 1.e-8 * (fabs ( Qcubed) + fabs( R * R)) )
      {
	double theta;
	if (Qcubed <= 0.0)theta = acos(1.0);
	else if (R / sqrt(Qcubed) > 1.0)theta = acos(1.0); 
	else if (R / sqrt(Qcubed) < -1.0)theta = acos(-1.0); 
	else theta = acos(R / sqrt(Qcubed));
	double sqrtQ = sqrt(Q);
	//	printf("sqrtQ = %12.5E teta=%12.5E a1=%12.5E\n",sqrt(Q),theta,a1);
	x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
	x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
	x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
      return (3);
      }

    /* Three real roots */
    if (d >= 0.0) {
      double theta = acos(R / sqrt(Qcubed));
      double sqrtQ = sqrt(Q);
      x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
      x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
      x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
      return (3);
    }
    
    /* One real root */
    else {
      printf("IMPOSSIBLE !!!\n");

      double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
      if (R > 0)
	e = -e;
      x[0] = (e + Q / e) - a1 / 3.;
      return (1);
    }
  }

#define MAXN 32
#define R(i,j)	result[n*(i)+(j)]
  
  long NullSpace(const double *a, double *result, double eps, long n)
  {
    int r[MAXN], c[MAXN];
    long i, j, k;
    int jj, kk, t;
    double max, temp;
    int ec;
    
    for (i = 0; i < n; i++)
      r[i] = c[i] = -1;			/* Reset row and column pivot indices */
    
    // copy the input matrix if not in place
    if (result != a) 
      for (i = 0; i < n*n; i++)  
	result[i] = a[i];
    // rest of algorithm is in place wrt result[]
    
    for (i = 0; i < n; i++) {
      /* Find the biggest element in the remaining submatrix
       * for the next full pivot.
       */
      max = 0.0;
      for (k = 0; k < n; k++) {
	if (r[k] < 0) {
	  for (j = 0; j < n; j++) {
	    if ((c[j] < 0) && ((temp = fabs(R(k, j))) > max)) {
	      kk = k;
	      jj = j;
	      max = temp;
	    }
	  }
	}
      }
      if (max < eps)
	break;		/* Consider this and all subsequent pivots to be zero */
      
      c[jj] = kk;					/* The row */
      r[kk] = jj;					/*	      and column of the next pivot */
      
      temp = 1.0 / R(kk, jj);
      R(kk, jj) = 1.0;
      for (j = 0; j < n; j++)		/* Should this be for j != jj ? */
	R(kk, j) *= temp;		/* Row equilibration */
      
      for (k = 0; k < n; k++) {	/* Row elimination */
	if (k == kk)
	  continue;			/* Don't do a thing to the pivot row */
	temp = R(k, jj);
	R(k, jj) = 0.0;
	for (j = 0; j < n; j++) {
	  R(k, j) -= temp * R(kk, j);	/* Subtract row kk from row k */
	  if (fabs(R(k, j)) < eps)
	    R(k, j) = 0.0;	/* Flush to zero if too small */
	}
      }
    }
    
    /* Sort into a truncated triangular matrix */
    for (j = 0; j < n; j++) {		/* For all columns... */
      while ((c[j] >= 0) && (j != c[j])) {
	for (k = 0; k < n; k++) {
	  if (r[k] < 0) {
	    /* Aha! a null column vector */
	    temp = R(k, j);	/* Get it on top */
	    R(k, j) = R(k, c[j]);
	    R(k, c[j]) = temp;
	  }
	}
	t = c[j];				/* Twiddle until pivots are on the diagonal */
	c[j] = c[t];
	c[t] = t;
      }
    }
    
    /* Copy the null space vectors into the top of the A matrix */
    ec = 0;
    for (k = 0; k < n; k++) {
      if (r[k] < 0) {
	R(k, k) = 1.0;			/* Set the pivot equal to 1 */
	if (ec != k) {
	  for (j = 0; j < n; j++) {
	    R(ec, j) = R(k, j);
	  }
	}
	ec++;
      }
    }
    /* The first  ec  rows of the matrix  a  are the vectors which are
     * orthogonal to the columns of the matrix  a.
     */
    return (ec);
  }    

  struct greater_abs
  {
    bool operator () (const double &a, const double &b)
    {
      return fabs(a) > fabs(b);
    }
  };

  long mTensor2::eigen (mVector e[3], double v[3]) const
  {            
    /// characteristic polynomial of T : find v root of
    /// v^3 - I1 v^2 + I2 T + I3 = 0
    /// I1 : first invariant , trace(T)
    /// I2 : second invariant , 1/2 (I1^2 -trace(T^2))
    /// I3 : third invariant , det T
    double I[4];
    I[3] = 1.0;
    I[2] = - trace();
    I[1] = 0.5 * (I[2]*I[2] - trace2());
    I[0] = - det();

    //    printf (" %lf x^3 +  %lf x^2 + %lf x + %lf = 0\n",
    //    	  I[3],I[2],I[1],I[0]);

    long nbEigen = FindCubicRoots (I,v);

    std::sort(v,v+3, greater_abs() );
    
    //    printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);

    double result[12];
    int nb_vec=0;
    while(1)
      {
	double a[9] = {pos[0][0]-v[nb_vec],pos[0][1],pos[0][2],
		       pos[1][0],pos[1][1]-v[nb_vec],pos[1][2],
		       pos[2][0],pos[2][1],pos[2][2]-v[nb_vec]};
	
	double eps = 1.e-3;
	int nb = 0;
	while (1)
	  {
	    nb = NullSpace (a,result,eps,3);
	    if (nb != 0)break;
	    eps *= 2.0;
	  }
	int kk=0;
	for (int i=nb_vec;i<nb+nb_vec;i++)
	  {
	    e[i] = mVector (result[0+kk*3],result[1+kk*3],result[2+kk*3]);
	    e[i].norm();
	    kk++;
	    if (i == 2)return nbEigen;
	  }
	nb_vec += nb;
	if (nb_vec == 3)return nbEigen;
	if (nb > 3)throw;
      }
    //    printf (" %lf x^3 +  %lf x^2 + %lf x + %22.15E = 0\n",
    //	  I[3],I[2],I[1],I[0]);
    throw;
  }


  long mTensor2::eigen2d (mVector e[3], double v[3]) const
  {            

    double a = 1.0;
    double b = - pos[0][0] - pos[1][1];
    double c = (pos[0][0] * pos[1][1] - pos[1][0] * pos[0][1]);

    e[2] = mVector(0,0,1);

    double delta = b*b - 4 * a * c;

    if (delta < 0)return 0;
    

    v[0] = (-b+sqrt(delta))/(2.*a);
    v[1] = (-b-sqrt(delta))/(2.*a);
    v[2] = 1.0;

    long nbEigen = 2;

    if (fabs(v[1]) > fabs(v[0]))
      {
	double temp = v[0];
	v[0] = v[1];
	v[1] = temp;
      }
    
    //    printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);

    double result[4];
    int nb_vec=0;
    while(1)
      {
	double a[4] = {pos[0][0]-v[nb_vec],pos[0][1],
		       pos[1][0],pos[1][1]-v[nb_vec]};	
	double eps = 1.e-8;
	int nb = 0;
	while (1)
	  {
	    nb = NullSpace (a,result,eps,2);
	    if (nb != 0)break;
	    eps *= 2.0;
	    //	    printf ("esp = %12.5E\n",eps);	
	  }
	int kk=0;
	for (int i=nb_vec;i<nb+nb_vec;i++)
	  {
	    e[i] = mVector (result[0+kk*2],result[1+kk*2],0.0);
	    e[i].norm();
	    kk++;
	  }
	nb_vec += nb;
	if (nb_vec == 2)return nbEigen;
	if (nb > 2)throw;
      }
    throw;
  }

/*
  mTensor2 mTensor2::transpose() const
  {
    mTensor2 tt;
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	{
	  tt(i,j) = pos[j][i];
	}
    return tt;
  }
*/

  void mTensor2::transpose()
  {
    for (int i=0;i<3;i++)
      for (int j=i+1;j<3;j++) {
                double temp = pos[j][i];
                pos[j][i] = pos[i][j];
                pos[i][j] = temp;
          }
  }

  mTensor2 mTensor2::operator !() const
  {
        mTensor2 B(*this);
        B.transpose();
        return B;
  }

  mTensor2 mTensor2::intersect2d (mTensor2 &t) const
  {
    
    mVector e1[3],e2[3];
    double v1[3],v2[3];

    eigen2d (e1,v1);
    t.eigen2d (e2,v2);
        
    if (fabs(v1[0]) > fabs(v2[0]))
      {
	mTensor2 TR (e1[0],e1[1],e1[2]);
	mTensor2 TRT = TR;
	TRT.transpose();
	mTensor2 D (0.0);
	D(0,0) = v1[0];
	D(2,2) = 1.0;
	double x = (t * e1[1]) * e1[1];
	if(fabs(x) > fabs(v1[1]))
	  D(1,1) = x;
	else
	  D(1,1) = v1[1];	
	return (TR*D)*TRT;
      }
    else
      {
	mTensor2 TR (e2[0],e2[1],e2[2]);
	mTensor2 TRT = TR;
	TRT.transpose();
	mTensor2 D (0.0);
	D(0,0) = v2[0];
	D(2,2) = 1.0;
	double x = ((*this) * e2[1]) * e2[1];
	if(fabs(x) > fabs(v2[1]))
	  D(1,1) = x;
	else
	  D(1,1) = v2[1];	
	return (TR*D)*TRT;
      }
    
  }

  mTensor2 mTensor2::operator * (const mTensor4 &other) const
    {
      mTensor2 m(0.0);
      for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	  m(i,j) = 0.0;
	  for(int k=0;k<3;k++){
	    for(int l=0;l<3;l++){
	      m(i,j) += pos[k][l] * other(i,j,k,l);
	    }
	  }
	}
      }
      return m;
    }


} // end of namespace

