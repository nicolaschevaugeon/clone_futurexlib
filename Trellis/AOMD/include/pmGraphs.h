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

#ifndef _MGRAPHS_
#define _MGRAPHS_

#include <iosfwd>


namespace AOMD {

class mMesh;

/**
  One graph contains vertices and edges, these are not vertices and edges of 
  the mMesh ! The idea is that the vertices of the graph are the mesh entities 
  of dimension from and the edges are entities of dimesion to.

  The graph is designed as a struct because it's a very internal structure.

  AOMD_graph::nn is the number of vertices of the graph
  
  AOMD_graph::xadj is a vector of "pointers", its dimension is (nn+1) 
  xadj[i] points to the first element of the adjacency of node i in adjncy. 
  The terminology is exactly the same as the one you can find in the METIS documentation.

  NOTE : EVERYTHING USES FORTRAN INDEXING so make j = xadj[k] - 1 to access adjncy[j]

  xadj [ 0 , ...,  i , ... , nn ]

  adjncy [ 0 , ... , xadj[i] - 1 , ... , xadj[i+1] - 1, ... ]
                     |__________________|
                            
                     indexes of vertices 
		       adjacent to i
   
  vwgt of size nn contains weights on vertices
  adjwgt of big size contains weights on edges

  In the constructor, you give a bool (doWeUseWeights) telling if you want to compute
  initial weights based on element sizes.
*/

struct AOMD_graph
{
  unsigned int tagDn;
  AOMD_graph(mMesh *, int from, int to , bool doWeUseWeights);
  AOMD_graph(mMesh *, int from, int to , bool doWeUseWeights, int level);
  AOMD_graph();
  virtual ~AOMD_graph();
  int nn;
  int *xadj;
  int *adjncy;
  int * vwgt, *adjwgt;
  void print(std::ostream &);
};

class AOMD_LoadBalancerCallbacks;
void computeNbNodesPerProc(mMesh *, int *nbNodesPerProc, int from, int level=-1);

/**
  Distributed graph, exactly the same as a mGraph except a new vector
  vtxdist that has to be the same in all processes.

  np is the number of processes
  vtxdist is of size np+1
  
  vertices are inexed, indexes have to be different on
  different processes, each process has a range of indexes.

  vtxdist[i] gives the initial vertex id on process i
*/


class AOMD_distributed_graph
{
 public:
  AOMD_distributed_graph(mMesh *, int from, int to , 
			 AOMD_LoadBalancerCallbacks *dm = nullptr,
			 int level = -1);
  AOMD_distributed_graph();
  virtual ~AOMD_distributed_graph();
  AOMD_graph *theGraph;
  int *vtxdist;
  void print(std::ostream &);
};



} // end of namespace

#endif
