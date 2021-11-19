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
#ifndef _AOMD_LOADBALANCER_H_
#define _AOMD_LOADBALANCER_H_

namespace AOMD {

class mEntity;
class AOMD_distributed_graph;
/**
   Load balancing callbacks for AOMD. User must
   derive a class form AOMD_LoadBalancerCallbacks
   if he wants to load balance a mesh. mMesh::loadBalance
   takes this callback as parameter.
*/
class AOMD_LoadBalancerCallbacks
{
public :
  /// from a given graph, we retrieve a partitionVector
  virtual void partition(AOMD_distributed_graph &theGraph , int *partitionVector) = 0;
  /// do we use a weighted graph
  virtual bool useWeights() const = 0;
  /// node weight for a given mEntity
  virtual float getWeight (mEntity *) const = 0;
  /// user can create data's attached to a given mesh entity
  /// data's have size "size" 
  virtual void * getUserData (mEntity *, int dest_proc, int &size) = 0; 
  /// user has to provide a way to delete its own data (delete or free).
  virtual void deleteUserData (void *) = 0; 
  /// user recieves its data's. mEntity is now the mesh entity on the 
  /// remote side.
  virtual void recieveUserData (mEntity *, int pid, int tag, void *buf) = 0;
  virtual int nbProcs() const = 0;
};

} // end of namespace

#endif




