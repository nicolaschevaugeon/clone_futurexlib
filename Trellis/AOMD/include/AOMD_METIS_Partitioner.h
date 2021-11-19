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
#ifndef _AOMD_METIS_PARTITIONER_H_
#define _AOMD_METIS_PARTITIONER_H_
#include "AOMD_LoadBalancer.h"
#ifdef _HAVE_PARMETIS_

namespace AOMD {
class AOMD_METISPartitioner : public AOMD_LoadBalancerCallbacks
{
  int NP;
public :
  /// Constructor takes the number of partitions as input
  AOMD_METISPartitioner (int);
  /// this function interfaces with Metis
  virtual void partition(AOMD_distributed_graph &theGraph , 
			 int *partitionVector);
  virtual bool useWeights() const {return false;}
  /// node weight for a given mEntity
  virtual float getWeight (mEntity *) const 
  { return 0.0;}
  /// user can create data's attached to a given mesh entity
  /// data's have size "size" 
  virtual void * getUserData (mEntity *, int dest_proc, int &size) 
  { return 0;} 
  /// user has to provide a way to delete its own data (delete or free).
  virtual void deleteUserData (void *) {} 
  /// user recieves its data's. mEntity is now the mesh entity on the 
  /// remote side.
  virtual void recieveUserData (mEntity *, int pid, int tag, void *buf) {} 
  virtual int nbProcs() const;
};

} // end of namespace 

#endif // PARALLEL
#endif //_AOMD_ZOLTAN_LOADBALANCER_H_
