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
#ifndef SIM
#ifndef H_AOMDLIST
#define H_AOMDLIST
#include <cstdio>
#include "AOMDfwd.h"
#include "AOMD.h"
#include "mEntity.h"
#include "GEntity.h"

/**
   This list is basically a vector of pointers
   with special constructor for adjacencies.
*/

namespace AOMD
{
  class fakeList
  {
    int ith;
    std::vector<void*> pointers;
  public :
    typedef std::vector<void*>::const_iterator iter;
    inline iter begin() const {return pointers.begin();}
    inline iter end()   const {return pointers.end()  ;}
    fakeList () : ith(0)
    {
    }
    fakeList (pEntity pe, int dim) : ith(0)
    {
      if(pe->getLevel() >= dim || pe->isAdjacencyCreated(dim))
	{
	  for(int i=0;i<pe->size(dim);i++) 
	    pointers.push_back((void*)pe->get(dim,i));
	}
      else
	{
	  mAdjacencyContainer upward;
	  pe->getHigherOrderUpward (dim,upward);
	  for(int i=0;i<upward.size();i++)
	    {
	      pointers.push_back(upward[i]);
	    }
	}
    }
    
    fakeList(GEntity* ge, int dim) : ith(0)
    {
      printf("This has not been implemented yet!\n");
      throw;
      /*      for(gEntity::iter it = ge->gadj_begin(); it != ge->gadj_end() ; ++it)
	{
	  if((*it)->getLevel() == dim)
	    pointers.push_back ((void*)*it);
	}
	*/
    }  
    inline int inlist (void *item)
    {
      std::vector<void*>::iterator it =
	std::find (pointers.begin(),pointers.end(),item);
      if(it == pointers.end()) return 0;
      return 1;
    }  
    inline void append (void *item)
    {
      pointers.push_back(item);
    }  
    inline void remove (void *item)
    {
      pointers.erase ( std::remove (pointers.begin(),
				    pointers.end(),item) , 
		       pointers.end () );
    }  
    inline void clear () 
    {
      pointers.clear();
    }
    inline int size () const 
    {
      return pointers.size();
    }

    inline void * operator [] (int i) 
    {
      return pointers[i];
    }

    inline void * next (void **restart)
    {
      
      static int totalgarbage;
      if(*restart == nullptr)
	{
	  ith=0;
	  *restart = &totalgarbage;
	}
      if(ith == size())
	{
	  ith = 0;
	  return nullptr;
	}
      return pointers[ith++];
    }  
    inline void * item (int i)
    {
      return pointers[i];
    }
  };
}
#endif
#endif
