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
#ifndef _MENTITY_CONTAINER_H_
#define _MENTITY_CONTAINER_H_
	
#include <vector>
#include <set>
#include <list>
#include <algorithm>

#ifndef SIM
#include "modeler.h"
#else
#include "SimModel.h"
#endif

const int _DIMS_ = 4;

#include "AOMDfwd.h"

#ifdef _NO_HASH_TABLE_
#include <set>
#else
#include <unordered_set>
#endif

namespace AOMD {
  /**
     Hash function for a mesh entity
  */
  class EntityHashKey 
    {
    public:
      int operator()(mEntity* ent2) const; 
    };

  /**
     Equal operator for 2 mesh entities
  */
  class EntityEqualKey 
    {
    public:
      bool operator()(mEntity* ent1, mEntity* ent2) const; 
    };

  /**
     Less than operator for 2 mesh entities
  */
  class EntityLessThanKey 
    {
    public:
      bool operator()(mEntity* ent1, mEntity* ent2) const;
    };

  /**
     A container for upward and downward entities.
     This is an internal of AOMD and that should not be
     accessed. User has only access to iterators
     through mEntity interface.

     New implementation, they are all vector<mEntity *>
  */

  class mAdjacencyContainer
    {
    private:
      /// all adjacency containers are sdt::vectors<mEntity*>
      /// Design has changed because std::set for upwards was
      /// costy, pretty slow and was not allowing inlining because
      /// of virtual functions
      std::vector<mEntity*> mContainer;
    public:
      /// iter is the const iterator that you should use to browse
      /// adjacency
      typedef std::vector<mEntity*>::const_iterator iter;
      /// Constructor, reserves given room
      inline mAdjacencyContainer(int i=1){mContainer.reserve(i);}
      /// nothing to be done here
      inline ~mAdjacencyContainer() = default;
      /// random access iterator
      inline mEntity* operator [] (int i); 
      /// add a mesh entity at the end of the container (push_back should be 
      /// more appropriate
      inline void add(mEntity* e);
      /// add a mesh entity to the container if it's not there yet
      inline void appendUnique(mEntity* e);
      /// deletes a mesh entity form the adjacency (do not delete the 
      /// mesh entity e itself !).
      inline void del(mEntity* e);
      /// same as operator []
      inline mEntity* get(int i) const;
      /// give how many adjacencies in the container
      inline int size() const;
      /// find if an adjacency is in the container
      inline mEntity* find(mEntity*) const;
      /// iterators
      inline iter begin() const {return mContainer.begin();}
      /// iterators
      inline iter end() const {return mContainer.end();}
      /// clear the container
      inline void clear();
    };

  inline mEntity*  mAdjacencyContainer::find(mEntity* ent) const 
    { 
      iter it = std::find (begin(),end(),ent);
      if(it == end())return nullptr;
      return *it;
    }

  inline void mAdjacencyContainer::clear() 
    {
      mContainer.clear();
    }
  
  inline mEntity* mAdjacencyContainer::operator[](int i) 
    { 
      return(mContainer[i]);
    }

  inline mEntity* mAdjacencyContainer::get(int i) const
    { 
      return(mContainer[i]);
    }

  inline void mAdjacencyContainer::del(mEntity* e)
    { 
      mContainer.erase ( std::remove (mContainer.begin(),mContainer.end(),e) , 
			 mContainer.end () );
    }

  inline void mAdjacencyContainer::add(mEntity* e)
    {
      mContainer.push_back(e);
    }

  inline void mAdjacencyContainer::appendUnique(mEntity* e)
    {
      if(!find(e))add(e);	
    }
  
  inline int mAdjacencyContainer::size() const 
    { 
      return mContainer.size();
    }

  /**
     A container for mesh entities (hash tables).
     This is an internal of AOMD and that should not be
     accessed. User has only access to iterators
     through mMesh interface.
  */

  class mMeshEntityContainer
    {
      public :
#ifdef _NO_HASH_TABLE_
      typedef std::set<mEntity*,EntityLessThanKey> CONTAINER;
#else
      typedef std::unordered_set<mEntity*,EntityHashKey,EntityEqualKey> CONTAINER;
#endif
      typedef CONTAINER::iterator iter;
    private:
      std::list<mIteratorNew*> currentIterators;
      CONTAINER *mEntities[_DIMS_];
    public:
      mMeshEntityContainer();
      virtual ~mMeshEntityContainer();
      iter begin(int what) const ;
      iter end(int what) const ;
      void add(mEntity* e);
      void del(mEntity* e);
      int  size(int what) const ;
      mEntity* find(mEntity*) const; 
      mIteratorNew* NewIterator (int what, int t=-2);
      mIteratorNew* NewIterator (int what, pGEntity);
      void DeleteIterator (mIteratorNew*);
    };

} // end of namespace
#endif

