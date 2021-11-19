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

#ifndef _MITERATORNEW_H_
#define _MITERATORNEW_H_

#include <list>
#include <vector>

#include "mEntityContainer.h"

class GEntity;

namespace AOMD
{
class mEntity;
class mMesh;

class mIteratorNew
{
  public:
#ifdef ITER_DEBUG
   static int id_generator;
   static std::vector<mIteratorNew *> all_iters;
#endif
   typedef std::forward_iterator_tag iterator_category;
   typedef class mEntity value_type;
   typedef class mEntity *pointer;
   typedef class mEntity &reference;
   typedef std::list<mEntity *>::iterator iter;
   typedef std::list<mEntity *>::const_iterator citer;

  private:
#ifdef ITER_DEBUG
   int id;
#endif
   int typ;  // type 1 : entityset iterator
             // type -2 ~ 0 : mesh iterator
   int dim;
   GEntity *ge;
   mMeshEntityContainer *allEntities;

   std::list<mEntity *> container;
   iter theIter;
// added for TSTT - seole Jul 09, 2004
#ifdef TSTT_
   int wSize;  // workset size
  public:
   // this one is added for TSTT entityset - seole 07/06/04
   // type 1 - entityset iterator
   mIteratorNew(std::list<mEntity *> &c, int type, int topo, int wsize = 1);
   // type 2 - mesh iterator
   mIteratorNew(mMesh *mesh, int type, int topo, int wsize = 1);

   int wrkSize() { return wSize; }
   int size() { return (int)container.size(); }
   int type() { return typ; }
#endif

  public:
   mIteratorNew(mMeshEntityContainer *c, int dim, int t);
   mIteratorNew(mMeshEntityContainer *c, int dim, GEntity *clas);
   ~mIteratorNew();

   mMeshEntityContainer *getEntities() { return allEntities; }
   friend class mMeshEntityContainer;
   void append(mEntity *e);
   void erase(mEntity *e);
   inline bool operator!=(mIteratorNew const &) const;
   inline bool operator==(mIteratorNew const &) const;
   inline mEntity *operator*() const;
   inline bool end() const;
   inline void next();
   inline void reset();
   inline mIteratorNew &operator++()
   {
      ++theIter;
      return *this;
   }
   inline mIteratorNew operator++(int)
   {
      mIteratorNew tmp = *this;
      ++(*this);
      return tmp;
   }
};

inline bool mIteratorNew ::end() const
{
   citer iend = theIter;
   return container.end() == iend;
}

inline void mIteratorNew ::next() { ++(*this); }

inline void mIteratorNew ::reset() { theIter = container.begin(); }

inline mEntity *mIteratorNew ::operator*() const { return *theIter; }

inline bool mIteratorNew ::operator!=(mIteratorNew const &other) const { return theIter != other.theIter; }

inline bool mIteratorNew ::operator==(mIteratorNew const &other) const { return theIter == other.theIter; }
}  // namespace AOMD
#endif
