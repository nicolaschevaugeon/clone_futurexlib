/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef _MESHINTERFACEXREGION2_
#define _MESHINTERFACEXREGION2_

#include <memory>

// AOMD
#include "mAOMD.h"
#include "mEdge.h"
#include "mFace.h"
#include "mRegion.h"
#include "mVertex.h"

// xfiles
#include "xRegion.h"

// xfastmarching
#include "linearalgebra3d.h"

namespace xfastmarching
{
template <class IT>
class citerobjecttovertexiteratorconvertor
{
  public:
   typedef const AOMD::mVertex *value_type;
   citerobjecttovertexiteratorconvertor(IT _it) : it(_it) {}
   citerobjecttovertexiteratorconvertor<IT> operator++(int)
   {
      citerobjecttovertexiteratorconvertor<IT> rit(it);
      ++(*this);
      return rit;
   }
   citerobjecttovertexiteratorconvertor<IT> &operator++()
   {
      ++it;
      return *this;
   }
   bool operator!=(const citerobjecttovertexiteratorconvertor &other) { return it != other.it; }
   void operator=(value_type v) { throw -1987; }
   const AOMD::mVertex *operator*() { return static_cast<const AOMD::mVertex *>(*it); }

  private:
   IT it;
};

template <class IT>
class entitytovertexiteratorconvertor
{
  public:
   typedef const AOMD::mVertex *value_type;
   entitytovertexiteratorconvertor(IT _it) : it(_it) {}
   entitytovertexiteratorconvertor<IT> operator++(int)
   {
      entitytovertexiteratorconvertor<IT> rit(it);
      ++(*this);
      return rit;
   }
   entitytovertexiteratorconvertor<IT> &operator++()
   {
      ++it;
      return *this;
   }
   bool operator!=(const entitytovertexiteratorconvertor &other) { return it != other.it; }
   void operator=(value_type v) { it = const_cast<AOMD::mVertex *>(v); }
   const AOMD::mVertex *operator*() { return static_cast<const AOMD::mVertex *>(*it); }

  private:
   IT it;
};

struct filterVertex
{
   bool operator()(const AOMD::mEntity *e) const
   {
      if (e->getLevel() < 1)
         return true;
      else
         return false;
   }
};

class meshinterfacexRegion
{
  public:
   meshinterfacexRegion(const xfem::xRegion &_reg) : reg(_reg), f_reg(nullptr) {}
   typedef AOMD::mVertex vertex;
   typedef AOMD::mEdge edge;
   typedef AOMD::mFace face;
   typedef AOMD::mRegion region;
   typedef AOMD::mEntity entity;
   typedef xfem::xMesh::partman_t partman_t;
   typedef xfem::keyManagerSendOrReceive key_manager_sor_t;
   typedef xfem::keyManagerSendOrReceive::information_key_t information_key_t;
   typedef xfem::xFilteredRegion<partman_t::c_iter_object_t, filterVertex> filtered_region_t;
   typedef filtered_region_t::FilterIter FilterIter;
   const xfem::xRegion &reg;
   std::unique_ptr<filtered_region_t> f_reg;
};

inline std::size_t getTag(const meshinterfacexRegion &mi, const std::string &tagname)
{
   return AOMD::AOMD_Util::Instance()->lookupMeshDataId(tagname.c_str());
};

inline std::size_t getNewTag(const meshinterfacexRegion &mi, const std::string &tagname)
{
   return AOMD::AOMD_Util::Instance()->newMeshDataId(tagname.c_str());
}

inline void releaseTag(const meshinterfacexRegion &mi, const std::string &tagname)
{
   int tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId(tagname.c_str());
   AOMD::AOMD_Util::Instance()->deleteMeshDataId(tag);
}

inline void releaseTag(const meshinterfacexRegion &mi, const std::size_t &tag)
{
   AOMD::AOMD_Util::Instance()->deleteMeshDataId(tag);
}

class AttachableDataPointer : public AOMD::mAttachableData
{
  public:
   AttachableDataPointer(void *pnt_) : pnt(pnt_) {}
   ~AttachableDataPointer() override = default;
   void *pnt;
};

inline void *getAttachedDataPointer(const meshinterfacexRegion &mi, AOMD::mEntity &en, size_t tag)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(en.getData(tag));
   if (dataatt)
      return dataatt->pnt;
   else
      return nullptr;
}

inline void deleteData(const meshinterfacexRegion &mi, AOMD::mEntity &en, size_t tag) { en.deleteData(tag); }
inline void attachDataPointer(const meshinterfacexRegion &mi, AOMD::mEntity &en, size_t tag, void *data)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(en.getData(tag));
   if (dataatt)
   {
      dataatt->pnt = data;
      return;
   }
   dataatt = new AttachableDataPointer(data);
   en.attachData(tag, dataatt);
}

inline void *getAttachedDataPointer(const meshinterfacexRegion &mi, AOMD::mVertex &v, size_t tag)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
      return dataatt->pnt;
   else
      return nullptr;
}

inline void attachDataPointer(const meshinterfacexRegion &mi, AOMD::mVertex &v, size_t tag, void *data)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
   {
      dataatt->pnt = data;
      return;
   }
   dataatt = new AttachableDataPointer(data);
   v.attachData(tag, dataatt);
}

inline void deleteData(const meshinterfacexRegion &mi, AOMD::mVertex &v, size_t tag) { v.deleteData(tag); }

inline void *getAttachedDataPointer(const meshinterfacexRegion &mi, AOMD::mFace &v, size_t tag)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
      return dataatt->pnt;
   else
      return nullptr;
}

inline void attachDataPointer(const meshinterfacexRegion &mi, AOMD::mFace &v, size_t tag, void *data)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
   {
      dataatt->pnt = data;
      return;
   }
   dataatt = new AttachableDataPointer(data);
   v.attachData(tag, dataatt);
}

inline void deleteData(const meshinterfacexRegion &mi, AOMD::mFace &v, size_t tag) { v.deleteData(tag); }

template <int DIM>
xtool::xRange<xfem::xIter> getEntityRange(const meshinterfacexRegion &mi)
{
   return xtool::make_range(mi.reg.begin(DIM), mi.reg.end(DIM));
}

template <class VERTEXINPUTITERATOR>
void getVerticesNeighbors(const meshinterfacexRegion &mi, const AOMD::mVertex &v, VERTEXINPUTITERATOR it)
{
   size_t nbedges = v.size(1);
   for (size_t i = 0; i < nbedges; ++i)
   {
      AOMD::mEntity *e = v.get(1, i);
      if (mi.reg.IsInRegion(e))
      {
         for (size_t k = 0; k < 2; ++k)
         {
            const AOMD::mVertex *vik = static_cast<const AOMD::mVertex *>(e->get(0, k));
            if (vik != &v) *it++ = vik;
         }
      }
   }
}

template <class EDGEINPUTITERATOR>
void getEdges(const meshinterfacexRegion &mi, const AOMD::mVertex &v, EDGEINPUTITERATOR it)
{
   size_t n = v.size(1);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(1, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mEdge *>(e);
         ++it;
      }
   }
}

template <class FACEINPUTITERATOR>
void getFaces(const meshinterfacexRegion &mi, const AOMD::mVertex &v, FACEINPUTITERATOR it)
{
   size_t n = v.size(2);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(2, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mFace *>(e);
         ++it;
      }
   }
}

template <class REGIONINPUTITERATOR>
void getRegions(const meshinterfacexRegion &mi, const AOMD::mVertex &v, REGIONINPUTITERATOR it)
{
   size_t n = v.size(3);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(3, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mRegion *>(e);
         ++it;
      }
   }
}

template <class FACEINPUTITERATOR>
void getFaces(const meshinterfacexRegion &mi, const AOMD::mEdge &v, FACEINPUTITERATOR it)
{
   size_t n = v.size(2);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(2, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mFace *>(e);
         ++it;
      }
   }
}

template <class REGIONINPUTITERATOR>
void getRegions(const meshinterfacexRegion &mi, const AOMD::mEdge &v, REGIONINPUTITERATOR it)
{
   size_t n = v.size(3);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(3, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mRegion *>(e);
         ++it;
      }
   }
}

template <class VERTEXINPUTITERATOR>
void getVertices(const meshinterfacexRegion &mi, const AOMD::mEdge &v, VERTEXINPUTITERATOR it)
{
   size_t n = v.size(0);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(0, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mVertex *>(e);
         ++it;
      }
   }
}

template <class VERTEXINPUTITERATOR>
void getVertices(const meshinterfacexRegion &mi, const AOMD::mFace &v, VERTEXINPUTITERATOR it)
{
   size_t n = v.size(0);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(0, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mVertex *>(e);
         ++it;
      }
   }
}

template <class VERTEXINPUTITERATOR>
void getVertices(const meshinterfacexRegion &mi, const AOMD::mRegion &v, VERTEXINPUTITERATOR it)
{
   size_t n = v.size(0);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(0, i);
      if (mi.reg.IsInRegion(e))
      {
         *it = static_cast<AOMD::mVertex *>(e);
         ++it;
      }
   }
}

inline double getX(const meshinterfacexRegion &mi, const AOMD::mVertex &v) { return v.point()(0); }

inline double getY(const meshinterfacexRegion &mi, const AOMD::mVertex &v) { return v.point()(1); }

inline double getZ(const meshinterfacexRegion &mi, const AOMD::mVertex &v) { return v.point()(2); }

inline int getId(const meshinterfacexRegion &mi, const AOMD::mVertex &v) { return v.getId(); }

template <class GEOMVECT>
void getCoord(const meshinterfacexRegion &mi, const AOMD::mVertex &v, GEOMVECT &x)
{
   std::cout << "getcoord not implemented for the used type " << std::endl;
   assert(0);
}

template <class T>
void getCoord(const meshinterfacexRegion &mi, const AOMD::mVertex &v, vector2d<T> &x)
{
   x = {v.point()(0), v.point()(1)};
}

template <class T>
void getCoord(const meshinterfacexRegion &mi, const AOMD::mVertex &v, vector3d<T> &x)
{
   x = {v.point()(0), v.point()(1), v.point()(2)};
}

inline int getNbVertex(const meshinterfacexRegion &mi) { return mi.reg.size(0); }

inline const AOMD::mVertex *getRemoteVertex(const meshinterfacexRegion &mi, const AOMD::mVertex &v, int remot_proc_id)
{
   auto po = mi.reg.getPartitionManager().getConstPartitionObject(static_cast<const AOMD::mEntity &>(v));
   return static_cast<const AOMD::mVertex *>(po.getRemoteObjectOn(remot_proc_id));
}

inline void setFilteredRegion(const meshinterfacexRegion &mi)
{
   const_cast<meshinterfacexRegion &>(mi).f_reg.reset(new meshinterfacexRegion::filtered_region_t(
       mi.reg.getPartitionManager().beginObject(), mi.reg.getPartitionManager().endObject(), filterVertex()));
   return;
}
inline citerobjecttovertexiteratorconvertor<meshinterfacexRegion::FilterIter> beginBnd(const meshinterfacexRegion &mi)
{
   return citerobjecttovertexiteratorconvertor<meshinterfacexRegion::FilterIter>(mi.f_reg->begin());
}
inline citerobjecttovertexiteratorconvertor<meshinterfacexRegion::FilterIter> endBnd(const meshinterfacexRegion &mi)
{
   return citerobjecttovertexiteratorconvertor<meshinterfacexRegion::FilterIter>(mi.f_reg->end());
}

}  // namespace xfastmarching
#endif
