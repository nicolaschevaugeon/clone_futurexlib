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
#include <cstdio>
#include <fstream>
#include "mAOMD.h"
#include "mMesh.h"
#include "mEntity.h"
#include "mBuildAdj.h"

namespace AOMD {

void PrintTopology (mMesh *m, std::ostream &o)
{
  printf("AOMD Warning: AOMD doesn't store model entity topology\n");
/*  mMesh::giter git = m->begin();
  mMesh::giter gitend = m->end();
  while(git != gitend)
  {
    gEntity *g = *git;
    o << "Ge(dim:" << g->getLevel() 
      << ",id:" << g->getId() 
      << ",tag:" << g->getTag() << ") <=> ";
    for (gEntity::iter it = g->gadj_begin();it != g->gadj_end();++it)
    {
      o << "Ge(dim:" << (*it)->getLevel() 
        << ",id" << (*it)->getId() 
        << ",tag:" << (*it)->getTag() << ") ";
    }
    o << "\n";
    git++;
  }
  */
}

  void buildForDimK ( mMesh *m, int k)
  {
    std::for_each(m->beginall(k),
		  m->endall(k),
		  createDownwardFunctor (k-1,0,m,false));
  }
  
  void connectForDimK (mMesh *m, int k)
  {
/*    printf("AOMD Warning: AOMD doesn't build model entity adjacency any more\n");
    mMesh::iterall it = m->beginall(k);
    mMesh::iterall itend = m->endall(k);
    while(it != itend)
    {
      mEntity *e = *it;
      GEntity *g = e->getClassification();
      if(e->isAdjacencyCreated(k-1))
      {
	int s = e->size(k-1);
	for(int i=0;i<s;i++)
	{
	  mEntity *sub = e->get(k-1,i);
	  GEntity *gsub = sub->getClassification();
	  //		if(abs(g->getLevel() - gsub->getLevel()) == 1)
	  {
	//		    g->gadj_add(gsub);
//		    gsub->gadj_add(g);
	  }
	}
      }  
      ++it;
    }
*/
  }
  
  void AOMD_Util::BuildGTopology ( mMesh *m )
  {
    int dim = m->getDim();
    {for(int i=2;i<=dim;i++) buildForDimK (m,i);}
    {for(int i=0;i< dim;i++) m->modifyState(i,i+1,true);}
    //printf("AOMD_Util::BuildTopology Warning: no more model entity adjacency from AOMD\n");
    //    {for(int i=1;i<=dim;i++) connectForDimK (m,i);}
    {for(int i=0;i< dim;i++) m->modifyState(i,i+1,false);}
    {for(int i=2;i<=dim;i++) m->modifyState(i,i-1,false);}
    m->reduceToMinumumRepresentation();
      
  }

  void AOMD_Util :: print_topology (const char *fname, mMesh *m)
  {
    std::ofstream ofs (fname);
    PrintTopology(m,ofs);
    ofs.close();
  }  

  void AOMD_Util :: BuildNullModel (mMesh *m)
  {
  }

}  // end of namespace

