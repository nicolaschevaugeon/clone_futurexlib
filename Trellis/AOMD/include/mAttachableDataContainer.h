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
#ifndef _MATTACHABLEDATACONTAINER_H_
#define _MATTACHABLEDATACONTAINER_H_

#include "mPoint.h"
#include "mTensor2.h"
#include "mVector.h"
#ifdef TSTT_
#include "mAOMD.h"
#endif
#include <algorithm>
#include <vector>

namespace AOMD
{
class mMesh;
class mMirrorVertex;
class mEntity;

/**
   Base class for attachable data's
*/
class mAttachableData
{
  public:
   virtual ~mAttachableData() = default;
};

/**
   Vector of ints as attachable data
*/
class mAttachableIntVector : public mAttachableData
{
  public:
   ~mAttachableIntVector() override = default;
   std::vector<int> v;
};

/**
   Vector of ints as attachable data
*/
class mAttachablePointVector : public mAttachableData
{
  public:
   ~mAttachablePointVector() override = default;
   std::vector<Trellis_Util::mPoint> v;
};

/**
   Int as attachable data
*/
class mAttachableInt : public mAttachableData
{
  public:
   ~mAttachableInt() override = default;
   int i;
};

/**
   Double as attachable data
*/
class mAttachableDouble : public mAttachableData
{
  public:
   ~mAttachableDouble() override = default;
   double d;
};

/**
   Point as attachable data
*/
class mAttachablePoint : public mAttachableData
{
  public:
   Trellis_Util::mPoint p;
};

/**
   Point as attachable data
*/
class mAttachable_mVector : public mAttachableData
{
  public:
   Trellis_Util::mVector v;
};

/**
   Point as attachable data
*/
class mAttachable_mTensor2 : public mAttachableData
{
  public:
   Trellis_Util::mTensor2 t;
};

/**
   Mesh as attachable data
*/
class mAttachableMesh : public mAttachableData
{
  public:
   mMesh *m;
};

/**
   Mesh as attachable data
*/
class mAttachableMirrorVertex : public mAttachableData
{
  public:
   mMirrorVertex *e;
};

/**
   mesh entity as attachable data
*/
class mAttachableEntity : public mAttachableData
{
  public:
   mEntity *e;
};

#ifdef TSTT_
class mAttachableOpaque : public mAttachableData
{
  public:
   char *o;
   int size;
};
#endif
/**
   Container for attachable data's Internal, mEntity and mMesh provides interfaces.
*/

class mAttachableDataContainer
{
  public:
   typedef std::pair<unsigned int, mAttachableData *> info;
   //      typedef std::map< unsigned int, mAttachableData *> container;
   typedef std::vector<info> container;
   typedef container::iterator iter_attachdata;
   typedef container::const_iterator citer_attachdata;

  private:
   container tab;

  public:
   mAttachableDataContainer() {}
   ~mAttachableDataContainer();
   inline void attachData(unsigned int, mAttachableData *);
   inline void deleteData(unsigned int);
   inline mAttachableData *getData(unsigned int) const;
   inline citer_attachdata begin_attachdata() const { return tab.begin(); };
   inline citer_attachdata end_attachdata() const { return tab.end(); };
   /// specific data types for mEntity
   inline void attachEntity(unsigned int, mEntity *);
   /// specific data types for int
   inline void attachInt(unsigned int, int);

   /// specific data types for int
   inline int getAttachedInt(unsigned int) const;
   inline int getAttachedInt(unsigned int, int *);

   /// specific data types for mEntity
   inline mEntity *getAttachedEntity(unsigned int);
   /// specific data types for double
   inline void attachDouble(unsigned int, double);

   /// specific data types for double
   inline double getAttachedDouble(unsigned int);
   inline int getAttachedDouble(unsigned int, double *);

   /// specific data types for vector
   inline void attachVector(unsigned int, const Trellis_Util::mVector &t);
   /// specific data types for vector
   inline Trellis_Util::mVector getAttachedVector(unsigned int);
   /// specific data types for tensor
   inline void attachTensor(unsigned int, const Trellis_Util::mTensor2 &t);
   /// specific data types for tensor
   inline Trellis_Util::mTensor2 getAttachedTensor(unsigned int);
#ifdef TSTT_
   inline void attachOpaque(unsigned int c, const char *o, int size);
   //  inline int getAttachedEntity(unsigned int c, mEntity* d);
   inline int getSizeAttachedOpaque(unsigned int c);
   inline int getAttachedOpaque(unsigned int c, char *d);
#endif
};

class equalInfoPred
{
   const unsigned int c;

  public:
   equalInfoPred(unsigned int i) : c(i) {}
   inline bool operator()(const mAttachableDataContainer::info &i) const { return (i.first == c); }
};

inline mAttachableDataContainer::~mAttachableDataContainer()
{
   if (!tab.size()) return;
   citer_attachdata it = begin_attachdata();
   citer_attachdata itEnd = end_attachdata();
   for (; it != itEnd; ++it)
   {
      mAttachableData *a = (*it).second;
      if (a)
      {
         delete a;
         a = nullptr;
      }
   }
}

inline mAttachableData *mAttachableDataContainer::getData(unsigned int c) const
{
   citer_attachdata it = std::find_if(tab.begin(), tab.end(), equalInfoPred(c));
   if (it == tab.end()) return nullptr;
   return (*it).second;
}

inline void mAttachableDataContainer::attachData(unsigned int c, mAttachableData *v) { tab.emplace_back(c, v); }

inline void mAttachableDataContainer::deleteData(unsigned int c)
{
#ifdef TSTT_
   if (AOMD_Util::Instance()->typeMeshDataId(c) == 3)
   {
      mAttachableOpaque *data = (mAttachableOpaque *)getData(c);
      delete[] data->o;
      delete data;
   }
   else
   {
#endif
      mAttachableData *data = getData(c);
      if (data) delete data;
#ifdef TSTT_
   }
#endif
   tab.erase(std::remove_if(tab.begin(), tab.end(), equalInfoPred(c)), tab.end());
}

inline void mAttachableDataContainer::attachInt(unsigned int c, int i)
{
   mAttachableInt *ai = (mAttachableInt *)getData(c);
   if (!ai)
   {
      ai = new mAttachableInt;
      attachData(c, ai);
   }
   ai->i = i;
}

inline int mAttachableDataContainer::getAttachedInt(unsigned int c) const
{
   mAttachableInt *ai = (mAttachableInt *)getData(c);
   if (!ai) return 0;
   return ai->i;
}

inline void mAttachableDataContainer::attachDouble(unsigned int c, double d)
{
   mAttachableDouble *ai = (mAttachableDouble *)getData(c);
   if (!ai)
   {
      ai = new mAttachableDouble;
      attachData(c, ai);
   }
   ai->d = d;
}

inline double mAttachableDataContainer::getAttachedDouble(unsigned int c)
{
   mAttachableDouble *ai = (mAttachableDouble *)getData(c);
   if (!ai) return 0;
   return ai->d;
}

inline void mAttachableDataContainer::attachEntity(unsigned int c, mEntity *e)
{
   mAttachableEntity *ai = (mAttachableEntity *)getData(c);
   if (!ai)
   {
      ai = new mAttachableEntity;
      attachData(c, ai);
   }
   ai->e = e;
}

inline mEntity *mAttachableDataContainer::getAttachedEntity(unsigned int c)
{
   mAttachableEntity *ai = (mAttachableEntity *)getData(c);
   if (!ai) return nullptr;
   return ai->e;
}

inline void mAttachableDataContainer::attachVector(unsigned int c, const Trellis_Util::mVector &d)
{
   mAttachable_mVector *ai = (mAttachable_mVector *)getData(c);
   if (!ai)
   {
      ai = new mAttachable_mVector;
      attachData(c, ai);
   }
   ai->v = d;
}

inline Trellis_Util::mVector mAttachableDataContainer::getAttachedVector(unsigned int c)
{
   mAttachable_mVector *ai = (mAttachable_mVector *)getData(c);
   if (!ai) return 0;
   return ai->v;
}

inline void mAttachableDataContainer::attachTensor(unsigned int c, const Trellis_Util::mTensor2 &d)
{
   mAttachable_mTensor2 *ai = (mAttachable_mTensor2 *)getData(c);
   if (!ai)
   {
      ai = new mAttachable_mTensor2;
      attachData(c, ai);
   }
   ai->t = d;
}

inline Trellis_Util::mTensor2 mAttachableDataContainer::getAttachedTensor(unsigned int c)
{
   mAttachable_mTensor2 *ai = (mAttachable_mTensor2 *)getData(c);
   if (!ai) return 0;
   return ai->t;
}

// added by E. Seol
// it returns 1 if succeed, 0 otherwise
inline int mAttachableDataContainer::getAttachedInt(unsigned int c, int *i)
{
   mAttachableInt *ai = (mAttachableInt *)getData(c);
   if (!ai) return 0;
   *i = ai->i;
   return 1;
}

// added by E. Seol
// it returns 1 if succeed, 0 otherwise
inline int mAttachableDataContainer::getAttachedDouble(unsigned int c, double *d)
{
   mAttachableDouble *ai = (mAttachableDouble *)getData(c);
   if (!ai) return 0;
   *d = ai->d;
   return 1;
}
#ifdef TSTT_
inline void mAttachableDataContainer::attachOpaque(unsigned int c, const char *o, int size)
{
   mAttachableOpaque *ai = (mAttachableOpaque *)getData(c);

   if (ai)
   {
      delete[] ai->o;
   }
   else
   {
      ai = new mAttachableOpaque;
      attachData(c, ai);
   }
   ai->size = size;
   ai->o = new char[size];
   memcpy(ai->o, o, size);
}

// this doesn't work...
/*  inline int mAttachableDataContainer::getAttachedEntity(unsigned int c, mEntity* e)
  {
    mAttachableEntity *ai = (mAttachableEntity *)getData(c);
    if(!ai)return 0;
    e = ai->e;
    return 1;
  }
*/

inline int mAttachableDataContainer::getSizeAttachedOpaque(unsigned int c)
{
   mAttachableOpaque *ai = (mAttachableOpaque *)getData(c);
   if (!ai) return 0;
   return ai->size;
}

inline int mAttachableDataContainer::getAttachedOpaque(unsigned int c, char *o)
{
   // we assume memory for char* o is allocated by application
   mAttachableOpaque *ai = (mAttachableOpaque *)getData(c);
   if (!ai) return 0;
   memcpy(o, ai->o, ai->size);
   return 1;
}
#endif

}  // namespace AOMD
#endif
