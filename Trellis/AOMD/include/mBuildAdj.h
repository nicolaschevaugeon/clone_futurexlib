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
#ifndef _MBUILDADJ_
#define _MBUILDADJ_

#include "mEntity.h"
#include "mMesh.h"

#ifndef SIM
#include "modeler.h"
#else
#include "SimModel.h"
#endif

namespace AOMD
{
  struct deleteAdjFunctor 
  {
    int dim;
    deleteAdjFunctor (int d) : dim (d) {}
    inline void operator () (mEntity *e) const
    {
      e->deleteAdjacencies(dim);
    }
  };

  struct createUpwardFunctor 
  {
    int dim;
    createUpwardFunctor (int d) : dim (d) {}
    inline void operator () (mEntity *ent) const
    {
      if (!ent->isAdjacencyCreated (dim))
	{
	  return;
	}
      int nb = ent->size(dim);
      for(int k = 0; k<nb ;k++)
	{
	  mEntity *sub = ent->get(dim,k);
	  std::list<mEntity *>leaves;
	  sub->getLeaves(leaves);
	  for(std::list<mEntity*>::const_iterator it2 = leaves.begin();
	     it2 != leaves.end();++it2)
	    {
	      mEntity *sup = *it2;
	      sup->appendUnique(ent);
	    }
	}
    }
  };

  struct createDownwardFunctor 
  {
    int i_dim, j_dim;
    mMesh *cont;
    bool force_create;
    bool create_upward_too;
    createDownwardFunctor (int i, int j, mMesh *m, bool f = true, bool up = false) 
      : i_dim (i), j_dim(j), cont(m), force_create(f),  create_upward_too(up) {}
    inline void operator () (mEntity *e) 
    {
      if((e)->isAdjacencyCreated(i_dim))return;
      for(int k=0;k<e->getNbTemplates(i_dim);k++)
	{
	  mEntity *t = e->getTemplate(k,i_dim,j_dim);
	  mEntity *q;
	  /// entity t already exists
	  if((q = (cont->find(t))))
	    {
	      if (GEN_type(e->getClassification()) < 
		  GEN_type(q->getClassification()) ) 
		{
		  q -> classify (e->getClassification());
		}
	      e->add(q);
	      if(create_upward_too)q->appendUnique(e);
	      if(q != t)delete t;
	    }
	  /// add the new one to the database
	  else if(force_create)
	    {
	      t->classify (e->getClassification());
	      //printf("yeeee\n");
	      e->add(t);
	      if(create_upward_too)t->appendUnique(e);
	      cont->add(t);
	    }
	  /// only add existing ones
	  else delete t;
	}
    }
  };
}// end of namespace AOMD
#endif
