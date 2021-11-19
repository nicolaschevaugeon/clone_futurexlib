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
 *  AOMD/cint/AOMD_Entity.cc
 *  Created by Seegyoung Seol, on Mon Dec 15 2003, 10:26:37 EDT
 *
 *  File Content: AOMD c-iterface for Entity
 *
 *************************************************************************** </i>*/
 
#ifndef SIM
#include "mAOMD.h"
#include "AOMD_cint.h"
#include "AOMD_Internals.h"
#include "AOMDfwd.h"
#include "AOMD_Defs.h"
#include "mEntity.h"
#include "mMesh.h"
#include "mPoint.h"
#include "mVertex.h"



#include <iostream>
#include <list>
#include <vector>

using namespace AOMD;
using std::cout;
using std::endl;
using std::vector;
using std::cerr;

bool EN_isGhost(pMesh mesh, pEntity ent)
{
  return false;

}

int EN_getWeight(pEntity pe, double* weight)
{
  return pe->getAttachedDouble(AOMD_Util::Instance()->getWeight(), weight);
}

void EN_setWeight(pEntity pe, double weight)
{
  pe->attachDouble(AOMD_Util::Instance()->getWeight(),weight);
}

std::string EN_getUidStr(pEntity ent)
{
  return ent->getUid();
}

int EN_owner(pEntity ent)
{
  return ent->getOwner();
}

void EN_print(pEntity ent)
{
  ent->print();
}


bool EN_duplicate(pEntity ent)
{
  bool ret = false;
  if (ent->getPClassification())
    ret = true;
//  cout<<"("<<ParUtil::Instance()->rank()<<") AOMD::EN_multiple("<<ent->getUid()<<")="<<ret<<endl;
  return ret;
}

/*************************************
  Remote Copy operators
**************************************/
pEntity EN_getCopy(pEntity ent,int pid)
{

  return (pEntity)nullptr;
 
}

void EN_getCopies(pEntity ent,
                        std::vector<std::pair<pEntity,int> >& remoteCopies)
{

}

void EN_addCopy(pEntity ent, int remotePid, pEntity remoteCopy)
{

}

void EN_clearCopies(pEntity ent)
{

}

/*************************************
  Partition Classification operators
**************************************/
AOMD::pmEntity* EN_getPClassification(pEntity ent)
{
  return ent->getPClassification();
}

void EN_setPClassification(pEntity ent, AOMD::pmEntity* pe)
{

}

void computeCrossProduct( double a[], double b[], double n[] ) 
{
  n[0] = a[1]*b[2] - b[1]*a[2] ;
  n[1] = a[2]*b[0] - a[0]*b[2] ;
  n[2] = a[0]*b[1] - a[1]*b[0] ;
}

double computeDotProduct( double a[], double b[] ) 
{
  double sum = 0.0 ;
  for( int i = 0 ; i < 3 ; ++i ) 
    sum += a[i]*b[i] ;
  return sum ;
}

double R_Volume(mEntity* tet) 
{
  double x[4];
  double y[4];
  double z[4];
  double xyz[3];
  for (int i=0; i<4; ++i)
  {
    V_coord((mVertex*)(tet->get(0,i)), xyz);
    x[i] =xyz[0];
    y[i] = xyz[1];
    z[i] = xyz[2];
  }
  double a[] = { x[1]-x[0], y[1]-y[0], z[1]-z[0] } ;
  double b[] = { x[2]-x[0], y[2]-y[0], z[2]-z[0] } ;
  double c[] = { x[3]-x[0], y[3]-y[0], z[3]-z[0] } ;
  double normal[3] ;
  computeCrossProduct( a, b, normal ) ;
  return computeDotProduct( normal, c )/6.0;
}

void R_center(pEntity ent, double xyz[3])
{
  Trellis_Util::mPoint pp;
  void *iter = nullptr;
  pEntity vt;
  int i;
  pPList vtxs = R_vertices((pRegion)ent,1);
  while (vt = (pEntity)PList_next(vtxs,&iter))
  {
    pp = ((mVertex*)vt)->point();
    xyz[0] += pp(0);
    xyz[1] += pp(1);
    xyz[2] += pp(2);
  }
  PList_delete(vtxs);
  if (ent->getType()==mEntity::TET)
  {
    for(i=0; i<3; i++) xyz[i] /= 4.0;
  }
  else if (ent->getType()==mEntity::PYRAMID)
  {
    for(i=0; i<3; i++) xyz[i] /= 5.0;
  }
  else // mEntity::PRISM or mEntity::HEX
  {
    for(i=0; i<3; i++) xyz[i] /= 6.0;
  }
}

void F_center(pEntity ent, double xyz[3])
{
  Trellis_Util::mPoint pp;
  void *iter = nullptr;
  pEntity vt;
  int i;
  pPList vtxs = F_vertices((pFace)ent,1);
  while (vt = (pEntity)PList_next(vtxs,&iter))
  {
    pp = ((mVertex*)vt)->point();
    xyz[0] += pp(0);
    xyz[1] += pp(1);
    xyz[2] += pp(2);
  }
  PList_delete(vtxs);
  if (ent->getType()==mEntity::TRI)
  {
    for(i=0; i<3; i++) xyz[i] /= 3.0;
  }
  else if (ent->getType()==mEntity::QUAD)
  {
    for(i=0; i<3; i++) xyz[i] /= 4.0;
  }
  else // mEntity::PRISM or mEntity::HEX
  {
   cerr<<"Unknown Face Type\n";
  }
}

void E_center(pEntity ent, double xyz[3])
{

  Trellis_Util::mPoint pp;
  pEntity vt;
  int i;
  for (i=0; i<2; ++i)
  {
    vt = ent->get(0,i);
    pp = ((mVertex*)vt)->point();
    xyz[0] += pp(0);
    xyz[1] += pp(1);
    xyz[2] += pp(2);
  }
   for(i=0; i<2; i++) xyz[i] /= 2.0;
}





#endif   /* ifndef SIM */
