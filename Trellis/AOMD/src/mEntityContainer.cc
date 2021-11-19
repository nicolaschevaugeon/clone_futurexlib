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
#include "mEntityContainer.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mException.h"
#include "mIteratorNew.h"

#include<cstdio>

using std::cout;
using std::list;

namespace AOMD {
  
  /************* mMeshEntityContainer *******************/
  mMeshEntityContainer::mMeshEntityContainer()
  {
    for(int i=0;i<_DIMS_;i++)mEntities[i] = new CONTAINER;
  }
  mMeshEntityContainer::~mMeshEntityContainer()
  {
    for(int i=0;i<_DIMS_;i++)delete mEntities[i];
  }
  
  mMeshEntityContainer::iter mMeshEntityContainer::begin(int what) const 
  {
#ifdef _DEBUG_
    if(what>3 || what < 0)
      throw new mException(__LINE__,__FILE__,"hey !! this is not a valid dimension for mesh entities iteration !!");
#endif
    return mEntities[what]->begin();
  }
  
  mMeshEntityContainer::iter mMeshEntityContainer::end(int what)   const 
  {
#ifdef _DEBUG_
    if(what>3 || what < 0)
      throw new mException(__LINE__,__FILE__,"hey !! this is not a valid dimension for mesh entities iteration !!");
#endif
    return mEntities[what]->end();
  }
  
  void mMeshEntityContainer::add(mEntity* e)
  { 
    mEntities[e->getLevel()]->insert(e);
    
    std::list<mIteratorNew*>::iterator it = currentIterators.begin();
    std::list<mIteratorNew*>::iterator itend = currentIterators.end();
    while (it != itend)
      {
	(*it)->append(e);
	++it;
      }
  }
  void mMeshEntityContainer::del(mEntity* e)
  { 
    mEntities[e->getLevel()]->erase(e);
    
    std::list<mIteratorNew*>::iterator it = currentIterators.begin();
    std::list<mIteratorNew*>::iterator itend = currentIterators.end();
    while (it != itend)
      {
	(*it)->erase(e);
	++it;
      }
  }
  int mMeshEntityContainer::size(int what) const 
  { 
    return mEntities[what]->size();
  } 
  
  mEntity* mMeshEntityContainer::find (mEntity *e) const 
  {
    iter it=mEntities[e->getLevel()]->find(e);
    if(it!=end(e->getLevel()))
      return(*it);
    else
      return(mEntity*)nullptr;
  }
  
  /*************KEYS **************/
  int EntityHashKey::operator()(mEntity* ent1) const 
  {
    return(ent1->getId()); 
  }
  
  bool EntityEqualKey::operator()(mEntity* ent1,mEntity* ent2) const 
  {
    return ent1->equal(ent2);
  }
  
  bool EntityLessThanKey::operator()(mEntity* ent1, mEntity* ent2) const 
  {
    return ent1->lessthan(ent2);
  }
  
  mIteratorNew *mMeshEntityContainer :: NewIterator (int what, int t)
  {
    mIteratorNew *it = new mIteratorNew(this,what,t);
    currentIterators.push_back(it);
    return it;
  }
  
  mIteratorNew *mMeshEntityContainer :: NewIterator (int what, GEntity *g)
  {
    mIteratorNew *it = new mIteratorNew(this,what,g);
    currentIterators.push_back(it);
    return it;
  }
  
  void mMeshEntityContainer :: DeleteIterator (mIteratorNew *it)
  {
    std::list<mIteratorNew*>::iterator iter = std::find(currentIterators.begin(),currentIterators.end(),it);
    delete *iter;
    currentIterators.erase(iter);
  }
} // end of namespace
