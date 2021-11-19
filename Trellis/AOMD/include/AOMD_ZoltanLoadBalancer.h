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
#ifndef _AOMD_ZOLTAN_LOADBALANCER_H_
#define _AOMD_ZOLTAN_LOADBALANCER_H_

#include "AOMD_LoadBalancer.h"
#include "ParUtil.h"

namespace AOMD {

/**
   A load balancer that uses Zoltan library.
   I have only implemented the interface for
   ParMetis and the Octree.
*/
class AOMD_ZoltanLoadBalancer : public AOMD_LoadBalancerCallbacks
{
public :
  /// Algorithms available
  typedef enum Algorithm {LDiffusion,GDiffusion,Remap,MLRemap,Random,Octree,Serial};
private:
  /// save the algorithm
  Algorithm theAlgorithm;
public :
  /// Constructor takes the algorithm as input
  AOMD_ZoltanLoadBalancer (Algorithm algo = Remap);
  /// this function interfaces with Zoltan
  virtual void partition(AOMD_distributed_graph &theGraph , 
			 int *partitionVector);
  /// change the algorithm
  void setAlgorithm (const Algorithm &algo);  
  virtual int nbProcs() const { return ParUtil::Instance()->size();}
};

} // end of namespace 


#endif //_AOMD_ZOLTAN_LOADBALANCER_H_

