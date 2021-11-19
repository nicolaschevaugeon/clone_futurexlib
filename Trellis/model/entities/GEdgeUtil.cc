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
#include "GEdgeUtil.h"
#include "SVector3.h"

using std::cerr;
double parAtDist2(const GEdge *e, double par, double dist, double tol, int dir);

double parAtDist(const GEdge *e, double par, double dist, double tol, int dir)
{
  AOMD::Range<double> prange = e->parBounds(0);
  double lastp, p=par;
  double sign, ds, pdist = 0;
  AOMD::SPoint3 pt1 = e->point(par);
  
  SVector3 der;
 
  sign = dir == 1 ? 1 : -1;
  int count = 0, error = 0;
  
  while( fabs(dist-pdist) > tol ){
    count++;
    // if we havent converged in 10 iterations, then something strange
    // is probably going on, so we'll bag it and do something that
    // has to work
    if(count > 10){
      error = 1;
      cerr << "parAtDist - didn't converge\n";
      break;
    }
    der = e->firstDer(par);
    ds = sign*(dist-pdist)/sqrt(der[0]*der[0] + der[1]*der[1] + der[2]*der[2]);
    lastp = p;
    p += ds;
    if(dir!=1){
      if(p < prange.low()){ // clamp to range of parameter space of edge
        p = prange.low();
	break;
      } else if( p > par)
	p = (lastp+par)/2.0;
    } else  { // dir == 1
      if (p > prange.high()) {
	p = prange.high();
	break;
      } else if (p < par)
	p = (lastp+par)/2.0;
    }
    AOMD::SPoint3 pt2 = e->point(p);
    pdist = sqrt( (pt1.x()-pt2.x())*(pt1.x()-pt2.x()) + 
		 (pt1.y()-pt2.y())*(pt1.y()-pt2.y()) +
		 (pt1.z()-pt2.z())*(pt1.z()-pt2.z()) );
  }
  if(error)
    p = parAtDist2(e,par,dist,tol,dir);
  return p;
}

double parAtDist2(const GEdge *e, double par, double dist, double tol, int dir)
{
  AOMD::Range<double> prange = e->parBounds(0);
  double p=par;
  double sign, ds, pdist = 0;
  AOMD::SPoint3 pt1 = e->point(par);

  SVector3 der;
 
  sign = dir == 1 ? 1 : -1;
  
  while( fabs(dist-pdist) > tol ){
    der = e->firstDer(par);
    // just step by the min of tol and distance to go
    double step = fabs(dist-pdist) > tol ? tol : dist-pdist;
    ds = sign*step/sqrt(der[0]*der[0] + der[1]*der[1] + der[2]*der[2]);
    p += ds;
    if(dir!=1){
      if(p < prange.low()){ // clamp to range of parameter space of edge
        p = prange.low();
	break;
      }
    } else  { // dir == 1
      if (p > prange.high()) {
	p = prange.high();
	break;
      }
    }
    AOMD::SPoint3 pt2 = e->point(p);
    pdist = sqrt( (pt1.x()-pt2.x())*(pt1.x()-pt2.x()) + 
		 (pt1.y()-pt2.y())*(pt1.y()-pt2.y()) +
		 (pt1.z()-pt2.z())*(pt1.z()-pt2.z()) );
  }
  return p;
}
