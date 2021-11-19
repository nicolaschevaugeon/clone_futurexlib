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
#include <iostream>
#include "pmGraphs.h"
#include "mMesh.h"
#include "mEntity.h"
#include "AOMD_LoadBalancer.h"
#include "ParUtil.h"
#include "mAOMD.h"



#include <cstdio>
#include <algorithm>
#include <map>

using namespace std;

namespace AOMD {

  /*
    Graph based on a certain level of refinement of the octree
    this is certaily better than the other way because the unbalance
    was much too large due to the strategy of keeping all subtrees
    on the same processor. Note that we will have to take into account 
    edge weights as well as vertex weights.
  */

  AOMD_graph::AOMD_graph (mMesh *m, int from, int to , bool w, int level)
  {
    int n = from;
    int nbTot = 0;    
    int numnod = 0;
  
    unsigned int tagId = AOMD_Util::Instance()->getId();
    unsigned int tagDn = AOMD_Util::Instance()->getDn();

    { // for win32 compiler
      for(mMesh::iterall it = m->beginall(n) ; it != m->endall(n) ; ++it)
	{
	  mEntity *e = *it;
	  // we only take this level
	  if(m->getRefinementLevel(e) == 0)
	    {	    
	      //	    nbTot++;
	      numnod++;
	      // access entities of dimension to (i.e. graph edges)
	      for(int i=0;i<e->size(to);i++)
		{
		  nbTot ++;
		}
	    }
	}
    }

    //  printf("sss proc %d %d %d\n",ParUtil::Instance()->rank(),numnod,nbTot);
    nn = numnod;

    /* allocation of the graph*/
    xadj   = new int[numnod+2];
    adjncy = new int[nbTot+1];    
    adjwgt = new int[nbTot+1];    

    /* temp vector */
    int icount = 0;
    int ecount = 0;  
    { 
      for(mMesh::iterall it = m->beginall(n) ; it!=m->endall(n) ; ++it)
	{
	  mEntity *e = *it;	  
	
	  if(m->getRefinementLevel(e) == 0)
	    {	    
	      std::set<int> adj;
	      std::map<int,int> wgt;
	      int id1 =  e->getAttachedInt(tagId);
	      for(int i=0;i<e->size(to);i++)
		{
		  mEntity *down = e->get(to,i);
		  list<mEntity*> leaves;
		  down->getLeaves(leaves);		

		  for(list<mEntity*>::const_iterator itl = leaves.begin();itl != leaves.end();++itl)
		    {
		      mEntity *leaf = *itl;
		      if(leaf->getData(tagDn))
			{
			  int idproc = leaf->getAttachedInt(tagDn);
			  //			  printf("inter proc found %d %d \n",idproc,id1);
			  adj.insert(idproc);
			  wgt[idproc] ++;
			}
		      for(mAdjacencyContainer::iter it2 = leaf->begin(n) ; it2!=leaf->end(n) ; ++it2)
			{
			  mEntity *up = *it2;
			  int id2 =  up->getAttachedInt(tagId) ;
			  if(id1 != id2)
			    {
			      adj.insert( id2 );
			      wgt[id2] ++;
			    }
			}
		    }
		}
	      xadj[ecount++] = icount+1;
	      for(std::set<int>::iterator its = adj.begin() ; its != adj.end() ; ++its) 
		{
		  adjncy[icount] = *its;
		  adjwgt[icount] = wgt[*its];
		  icount ++;
		}
	    }
	}
      xadj[ecount] = icount+1;
    }

    //  printf("proc %d nn %d icount %d ecount %d nbTot %d %d\n",ParUtil::Instance()->rank(),nn,icount,ecount,nbTot,numnod);

    vwgt=nullptr;
  }

  /** ------------------------------------------------------------------------
     This constructor is the one prensently used.     
   ---------------------------------------------------------------------------*/

  AOMD_graph::AOMD_graph (mMesh *m, int from, int to , bool w)
  {
    int n = from;

    unsigned int tagId = AOMD_Util::Instance()->getId();
    unsigned int tagDn = AOMD_Util::Instance()->getDn();

    int nbTot = 0;
    int numnod = 0;
    { // for win32 compiler

      // iterator on mesh entities of dimension n (i.e. graph vertices)
      for(mMesh::iter it = m->begin(n) ; it!=m->end(n) ; ++it)
	{
	  // access what is in the iterator
	  mEntity *e = *it;
	  numnod++;
	  // access entities of dimension to (i.e. graph edges)
	  for(int i=0;i<e->size(to);i++)
	    {
	      // down is ith the entity
	      mEntity *down = e->get(to,i);
	      list<mEntity*> leaves;
	      // it can be an octree-like structure, we get all the leaves !!!
	      down->getLeaves(leaves);
	      for(list<mEntity*>::const_iterator itl = leaves.begin();itl != leaves.end();++itl)
		{
		  mEntity *leaf = *itl;
		  if(leaf->getData(tagDn))nbTot++;
		  for(mAdjacencyContainer::iter it2 = leaf->begin(n) ; it2!=leaf->end(n) ; ++it2)
		    {
		      mEntity *up = *it2;
		      if(up != e)
			{
			  nbTot ++;
			}
		    }
		}
	    }
	}
    }
    
    nn = numnod;
    
    xadj   = new int[numnod+2];
    adjncy = new int[nbTot+1];    
    adjwgt = new int[nbTot+1];    
    
    int icount = 0;
    int ecount = 0;
    
    { 
      for(mMesh::iter it = m->begin(n) ; it!=m->end(n) ; ++it)
	{
	  mEntity *e = *it;
	  std::set<int> adj;
	  std::map<int,int> wgt;
	  int id1 =  e->getAttachedInt(tagId);
	  //      printf("node %d ",e->getAttachedInt("id"));e->print();
	  for(int i=0;i<e->size(to);i++)
	    {
	      mEntity *down = e->get(to,i);
	      list<mEntity*> leaves;
	      down->getLeaves(leaves);

	      for(list<mEntity*>::const_iterator itl = leaves.begin();itl != leaves.end();++itl)
		{
		  // for parallel purpose, we add the id of distributed  neighbor here
		  mEntity *leaf = *itl;

		  if(leaf->getData(tagDn))
		    {
		      int idproc = leaf->getAttachedInt(tagDn);
		      adj.insert(idproc);
		      wgt[idproc] = 1;
		      //		      printf("inter proc found %d %d \n",idproc,id1);
		    }

		  for(mAdjacencyContainer::iter it2 = leaf->begin(n) ; it2!=leaf->end(n) ; ++it2)
		    {
		      mEntity *up = *it2;
		      int id2 =  up->getAttachedInt(tagId) ;
		      if(id1 != id2)
			{
			  adj.insert( id2 );
			  if(up->root() == e->root())
			    wgt[id2]=10000;
			  else
			    wgt[id2]=1;
			}
		    }
		}
	    }
	  xadj[ecount++] = icount+1;      
	  for(std::set<int>::iterator its = adj.begin() ; its != adj.end() ; ++its) 
	    {
	      adjncy[icount] = *its;
	      adjwgt[icount] = wgt[*its];
	      icount ++;
	    }
	}
    }
    xadj[ecount] = icount+1;

    //  printf("proc %d nn %d icount %d ecount %d nbTot %d %d\n",ParUtil::Instance()->rank(),nn,icount,ecount,nbTot,numnod);

    //  print(std::cout);
    vwgt=nullptr;
  }

  AOMD_graph :: ~AOMD_graph()
  {
    delete [] xadj;
    delete [] adjncy;
    if(vwgt)delete [] vwgt;
    if(adjwgt)delete [] adjwgt;
  }

  void AOMD_graph::print( std::ostream & o)
  {
    o << "graph with " << nn << " nodes\n";

    o << "xadj ";
    for(int i=0;i<nn+1;i++)o << xadj[i] << " ";
    o << "\n";
    o << "adjncy ";
    for(int j=0; j<xadj[nn]-1;j++)
      {
	o << adjncy[j] << " ";
      }
    o << "\n";
    if(vwgt)
      {
	o << "vwgt ";
	for(int j=0; j<nn;j++)
	  {
	    o << vwgt[j] << " ";
	  }
	o << "\n";
      }
    if(adjwgt)
      {
	o << "adjwgt ";
	for(int j=0; j<xadj[nn]-1;j++)
	  {
	    o << adjwgt[j] << " ";
	  }
	o << "\n";
      }
  }



  AOMD_distributed_graph :: AOMD_distributed_graph(mMesh *m, int from, int to,
						   AOMD_LoadBalancerCallbacks *dm,
						   int level)  
  {
    /** Hypothesis is that nodes of the graph are mesh
	entities of dimension from and edges are of dimension
	to */
    int n = from;
    vtxdist = new int [ParUtil::Instance()->size()+1];
    int *nbNodesPerProc = new int [ParUtil::Instance()->size()];

    unsigned int tagId = AOMD_Util::Instance()->getId();
    unsigned int tagDn = AOMD_Util::Instance()->getDn();

    //printf("computing distributed graph (p%d)\n",ParUtil::Instance()->rank());

    /// We compute how many nodes of the graph per processor
    ParUtil::Instance()->Barrier(__LINE__,__FILE__);
    computeNbNodesPerProc(m,nbNodesPerProc,from,level);

    /// We then compute a vector that tells the range
    /// of nodes on processor i : 
    /// range of proc i is [vtxdist[i] ... vtxdist[i+1]]
    int init = 1;
    vtxdist[0] = 1;
    int total = 1;
    for (int i=0;i<ParUtil::Instance()->size();i++)
      {
	total += nbNodesPerProc[i];
	vtxdist[i+1] = total;      
      }
    
    delete [] nbNodesPerProc;
  

    /// Print vtxdist on the stdout
    /*if(ParUtil::Instance()->master())
      {
	ParUtil::Instance()->Msg(ParUtil::WARNING,"vtxdist(%d) = (",ParUtil::Instance()->rank());
	for (int i=0;i<=ParUtil::Instance()->size();i++)
	  {
	    if(i)
	     ParUtil::Instance()->Msg(ParUtil::WARNING,",%d",vtxdist[i]);
	    else
	      ParUtil::Instance()->Msg(ParUtil::WARNING,"%d",vtxdist[i]);
	  }
	ParUtil::Instance()->Msg(ParUtil::WARNING,")\n");	   
      }
    */
    /// We will tag all nodes with id's in this range
    init =vtxdist[ParUtil::Instance()->rank()];  
    /**
       if level < 0 then the graph is composed of all leaves
       of the mesh.

       else the graph is composed of elements of refinement
       level level. level=0 means all roots !!
    */

    if(level < 0)
      {
	for(mMesh::iter it = m->begin(n) ; it!=m->end(n) ; ++it)
	  {	  
	    (*it)->attachInt(tagId,init++);
	    if ((*it)->getAttachedInt(tagId) != init-1)throw;
	  }
      }
    else
      {
	for(mMesh::iterall it = m->beginall(n) ; it!=m->endall(n) ; ++it)
	  {
	    mEntity *e = *it;
	    // we only take this level	  
	    if(m->getRefinementLevel(e) == level)
	      {	    
		int id = init++;
		e->attachInt(tagId,id);
		list<mEntity*> leaves;
		e->getLeaves(leaves);		
		for(list<mEntity*>::const_iterator itl = leaves.begin();itl != leaves.end();++itl)
		  {
		    mEntity *leaf = *itl;
		    leaf->attachInt(tagId,id);
		  }
	      }
	  }
      }


    /// add neighbors through inter processor boundaries  

    
    // create the serial graph with inter processor links !!!
    if(level < 0)
      theGraph = new AOMD_graph (m,from,to,false);    
    else
      theGraph = new AOMD_graph (m,from,to,false,level);    
  


    /// add vertex weights to the graph
    {
      int initial_id = vtxdist[ParUtil::Instance()->rank()];
      if(theGraph->vwgt)delete [] theGraph->vwgt; 
      theGraph->vwgt = new int [theGraph->nn];
      for(int i=0;i<theGraph->nn;i++)theGraph->vwgt[i] = 1;
      for(mMesh::iter it = m->begin(n) ; it!=m->end(n) ; ++it)
	{
	  int id = (*it)->getAttachedInt(tagId) - initial_id;	
	  if(dm->useWeights())
	    theGraph->vwgt[id] *= (int) (dm->getWeight((*it)));
	}
      int sumw = 0;
      {      for(int i=0;i<theGraph->nn;i++)sumw += theGraph->vwgt[i];}
      //      ParUtil::Instance()->Msg(ParUtil::WARNING,"prw = %d\n",sumw);
    }
    


    //printf(" (%d) mGraph done...\n",ParUtil::Instance()->rank());
    for(mMesh::iter it = m->begin(n-1) ; it!=m->end(n-1) ; ++it)(*it)->deleteData(tagDn);
  }


  void AOMD_distributed_graph :: print(ostream &o)
  {
    printf("on processor %d :\n",ParUtil::Instance()->rank());
    theGraph->print(o);
    o << "vtxdist : ";
    for(int i=0;i<=ParUtil::Instance()->size();i++)o << vtxdist[i] << " ";
    o << "\n";
  }



  AOMD_distributed_graph :: ~AOMD_distributed_graph()
  {
    delete theGraph;
    delete [] vtxdist;
  }

  void computeNbNodesPerProc(mMesh *m, int *nbNodesPerProc, int n, int level)
  {
    int *sendcounts = new int[ParUtil::Instance()->size()];
    for(int i=0;i<ParUtil::Instance()->size();i++)sendcounts[i]=0;
  
    int numnod = 0;
    if(level < 0)
      {
	for(mMesh::iter it = m->begin(n) ; it!=m->end(n) ; ++it)
	  {
	    numnod ++;      	
	  }
      }
    else
      {
	for(mMesh::iterall it = m->beginall(n) ; it!=m->endall(n) ; ++it)
	  {
	    if(m->getRefinementLevel(*it) == 0)
	      numnod ++;      	
	  }
      }

    nbNodesPerProc[0] = numnod;
    delete[] sendcounts;

  }
} // end of namespace
