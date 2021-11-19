/******************************************************************************** 

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written 
   and maintained by the Scientific Computation Research Center (SCOREC) at 
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*********************************************************************************/

/*<i> ****************************************************************************
 *
 *  AOMD/cint/PARALLEL.cc
 *  Created by Seegyoung Seol, on Mon Dec 15 2003, 10:26:37 EDT
 *
 *  File Content: AOMD c-iterface for PARALLEL
 *
 *************************************************************************** </i>*/
 
#ifndef SIM
#include <algorithm>

#include "AOMD_cint.h"
#include "ParUtil.h"



using namespace AOMD;
using std::cout;
using std::endl;
using std::vector;

double P_getMaxDbl(double num)
{
  if (P_size()==1) return num;
  throw;

}

double P_getMinDbl(double num)
{
  if (P_size()==1) return num;
  throw;
}

double P_getSumDbl(double num)
{
  if (P_size()==1) return num;
  throw;
}

double P_wtime()
{return AOMD::ParUtil::Instance()->wTime(); }

void P_unifyMinMax(double* max, double* min)
{

}

void P_mergeArray(std::vector<int>& vec)
{

} 

int P_getMaxInt(int num)
{
  if (P_size()==1)    return num;
  throw;
}

int P_getMinInt(int num)
{
  if (P_size()==1) return num;
  throw;
}

int P_getSumInt(int num)
{
  if (P_size()==1) return num;
  throw;
}


void P_barrier()           //synchronization
{

}

int P_size()
{

  return 1;

}
 
int P_pid()
{

  return 0;

}

void P_getGlobal(int in, vector<int>& output)
{

  output.push_back(in);

}
 
#endif   /* ifndef SIM */
