/**************************************************************************** 

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

*****************************************************************************/
#ifndef _AOMD_OWNERMANAGER_H_
#define _AOMD_OWNERMANAGER_H_

/*
A container class for managing remote copies of
mesh entities in different partition boundaries
*/

#include <map>
#include <vector>
#include "AOMD_SharedInfo.h"
#include "mMesh.h"

namespace AOMD {

class mEntity;

class AOMD_OwnerManager 
{
  std::multimap <int,int> *theTopology;
  std::multimap <mEntity*,AOMD_SharedInfo> *allSharedInfos;
 public :
  // Added by A.S. 22/11/10 to suply the same service with mesh and submeseh in xRegion method xPartitionBoundaryInfo
  const  std::multimap <mEntity*,AOMD_SharedInfo> & getAllSharedInfo(){return * allSharedInfos;}
  // end addon
  typedef std::multimap<mEntity*,AOMD_SharedInfo>::const_iterator iter; 
  typedef std::multimap<int,int>::const_iterator iter_top; 
  AOMD_OwnerManager ();
  ~AOMD_OwnerManager ();
  void addRemoteCopy(mEntity* here, int pid, mEntity *there, bool p);
  // delete all the remote copies of this entity
  void deleteRemoteCopy(mEntity*);
  std::list<AOMD_SharedInfo> getSharedInfo(mEntity * e);
  void clearRemoteCopies();
  void addNeighbor(int, int otherside);
  bool isPeriodic(mEntity *);
  bool amITheOwnerOf (mEntity *);
  iter_top begin_top(int);
  iter_top end_top(int);
  iter_top begin_top();
  iter_top end_top();
  iter begin(mEntity*e);
  iter end(mEntity*e);
  iter begin();
  iter end();
  int count(mEntity*e);
  int count(int);
  void print();
  int find ( int nbProcs, int *procs);
  int create ( int nbProcs, int *procs);
  void getPartitionIds ( std::vector<int> &ids );
  void setPoorPid_asOwner(mMesh* mesh);
};

} // end of namespace

#endif

