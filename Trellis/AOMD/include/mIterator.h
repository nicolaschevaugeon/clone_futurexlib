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
#ifndef _MITERATOR_H_
#define _MITERATOR_H_

#include "mEntityContainer.h"
#include "mEntity.h"
#include <iterator>
#include <iostream>

#ifndef SIM
#include "modeler.h"
#else
#include "SimModel.h"
#endif

namespace AOMD {

typedef mEntity* EntityPtr_;

class mIterator: public std::iterator<std::forward_iterator_tag,EntityPtr_, int,const EntityPtr_ *,const EntityPtr_ &> {

public:
    typedef mMeshEntityContainer::iter iter;
  protected :
    iter theBeg;
    iter theEnd;
    iter theIter;
  public:
    inline mIterator(const iter &begin, const iter &end);
    virtual ~mIterator();
    inline bool operator != (mIterator const &) const;
    inline bool operator == (mIterator const &) const;
    inline reference operator *() const;
    inline bool end() const {return theIter == theEnd;}
    // This is being used for the C-interface. Don't use
    // it in C++! 
    virtual void next() = 0;
    inline void reset() {theIter = theBeg;}
  };

  class mFullIterator : public mIterator 
  {
  public:
    inline mFullIterator(const iter &begin, const iter &end);
    inline mFullIterator& operator++()
    {
      ++theIter;
      return *this;
    }
    inline mFullIterator operator++(int)
    {
      mFullIterator tmp = *this;
      ++(*this);
      return tmp;
    }
    void next() override;
  };

  class mLeavesIterator : public mIterator 
  {
  public:
    inline mLeavesIterator(const iter &begin, const iter &end);
    inline mLeavesIterator& operator++()
    {
      while(1)
	{
	  ++theIter;
	  if(theIter == theEnd)break;
	  mEntity *e = *theIter;
	  if(!e->isAdjacencyCreated(e->getLevel()))break;
	}
      return *this;
    }
    inline mLeavesIterator operator++(int)
    {
      mLeavesIterator tmp = *this;
      ++(*this);
      return tmp;
    }
    void next() override;
  };

  class mClassIterator : public mIterator 
  {
    int iClass;
    int iOnWhat;
  public:
    inline mClassIterator(const iter &begin, const iter &end, int i, int w);
    inline mClassIterator& operator++()
    {
      while(1)
	{
	  ++theIter;
	  if(theIter == theEnd)break;
	  mEntity *e = *theIter; 
	  pGEntity g = e->getClassification();
	  if(!e->isAdjacencyCreated(e->getLevel()))
	    if(GEN_tag(g) == iClass && GEN_type(g) == iOnWhat)
	      break;
	}
      return *this;
    }
    inline mClassIterator operator++(int)
    {
      mClassIterator tmp = *this;
      ++(*this);
      return tmp;
    }
    void next() override;
  };

  class mClassClosureIterator : public mIterator 
  {
    int iClass;
    int iOnWhat;
    mMeshEntityContainer *cont;
  public:
    inline mClassClosureIterator(const iter &begin, const iter &end, int i, int w);
    inline mClassClosureIterator& operator++();
    void next() override;
  };

  inline mIterator :: mIterator (const iter &begin, const iter &end)
    : theBeg(begin), theEnd (end), theIter(begin) 
  {
  }

  inline mFullIterator :: mFullIterator (const iter &begin, const iter &end)
    : mIterator(begin,end)
  {
  }

  inline mLeavesIterator :: mLeavesIterator (const iter &begin, const iter &end)
    : mIterator(begin,end)
  {
    while(1)
      {
	if(theIter == theEnd)break;
	mEntity *e = *theIter;
	if(!e->isAdjacencyCreated(e->getLevel()))break;
	++theIter;
      }
  }

  inline mClassIterator :: mClassIterator (const iter &begin, const iter &end, int c, int w)
    : mIterator(begin,end), iClass(c),iOnWhat(w)
  {
    while(1)
      {
	if(theIter == theEnd)break;
	mEntity *e = *theIter;
	pGEntity g = e->getClassification();
	if(!e->isAdjacencyCreated(e->getLevel()))
	  if(GEN_tag(g) == iClass && GEN_type(g) == iOnWhat)break;
	++theIter;
      }
  }

  inline mClassClosureIterator :: mClassClosureIterator (const iter &begin, const iter &end, int c, int w)
    : mIterator(begin,end), iClass(c),iOnWhat(w)
  {
    throw 1;
  }

  inline mIterator::reference mIterator :: operator * () const
  {
    return *theIter;
  }

  inline bool mIterator :: operator != (mIterator const & other) const
  {
    return theIter != other.theIter;
  } 

  inline bool mIterator :: operator == (mIterator const & other) const
  {
    return theIter == other.theIter;
  } 
} // end of namespace
#endif
