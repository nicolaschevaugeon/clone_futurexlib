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
#include <algorithm>
#include <cstdio>
#include <list>
#include "mIteratorNew.h"
#include "mEntity.h"
#include "mMesh.h"
#include <cassert>

using std::list;
using std::cout;
using std::endl;

namespace AOMD {

#ifdef ITER_DEBUG
int mIteratorNew::id_generator = 0;
vector<mIteratorNew*> mIteratorNew::all_iters;
#endif

mIteratorNew::~mIteratorNew ()
{
#ifdef ITER_DEBUG
  cout<<"-- mIterator "<<id<<" is deleted\n";
#endif
}
  
mIteratorNew::mIteratorNew (mMeshEntityContainer *c,
			      int what,
			      int t)
    : allEntities(c),typ(t),dim(what),ge(nullptr)
#ifdef TSTT_
    , wSize(1)
#endif
{
#ifdef ITER_DEBUG
    id = id_generator++;
    all_iters.push_back(this);
    cout<<"-- mIterator "<<id<<" is created\n";
#endif

    mMeshEntityContainer::iter b = c->begin(what);
    mMeshEntityContainer::iter e = c->end(what);
    switch (t)
      {
	/// all
      case -2:
{	for(mMeshEntityContainer::iter it = b;it != e;++it)
	  container.push_back(*it);}
	break;
	/// leaves
      case -1:
{	for(mMeshEntityContainer::iter it = b;it != e;++it)
	  if(!(*it)->isAdjacencyCreated((*it)->getLevel()))
	    container.push_back(*it);}
	break;
	/// roots
      case 0:
{	for(mMeshEntityContainer::iter it = b;it != e;++it)
	  if(!(*it)->parent())container.push_back(*it);}
	break;
      default :
	throw;
      }
    //    printf ("an iterator of type %d (dim %d) has been created, temp.size = %d\n",t,dim,container.size());
    theIter = container.begin();
  }

  mIteratorNew::mIteratorNew (mMeshEntityContainer *c,
			      int what,
			      GEntity *g)
    : allEntities(c),typ(0),dim(what),ge(g)
#ifdef TSTT_
    , wSize(1)
#endif
  {
#ifdef ITER_DEBUG
    id = id_generator++;
    all_iters.push_back(this);
    cout<<"-- mIterator "<<id<<" is created\n";
#endif
    
    mMeshEntityContainer::iter b = c->begin(what);
    mMeshEntityContainer::iter e = c->end(what);
    for(mMeshEntityContainer::iter it = b;it != e;++it)
      if((*it)->getClassification() == g)
        container.push_back(*it);
    theIter = container.begin();
  }

  void mIteratorNew::append (mEntity *e)
  {
    int d = e->getLevel();
    if (dim == d && typ == -2)
      {
	//	printf ("adding to the current iterator (dim = %d)\n",dim);
	container.push_back(e);
      }
  }

void mIteratorNew::erase (mEntity *e)
{
  int d = e->getLevel();
  if (dim == d && typ == -2)
  {
    //	printf ("i'm an iterator of dimension %d\n",dim);
    if (*theIter == e)
    {
      iter oldIter = theIter;
      ++theIter;
      container.erase(oldIter);
    }
    else
    {
      iter oldIter = theIter;
      oldIter --;
      if (*oldIter == e)
        container.erase(oldIter);
      else
      {
 	iter it = std::find(container.begin(), container.end(),e);
	if (it != container.end())
	{
	  cout<<"AOMD WARNING: "<<(*it)->getUid()<<" deleted from mIterator\n";
	  container.erase(it);
	}
      }
    }
  }
}
  
} // end of namespace
