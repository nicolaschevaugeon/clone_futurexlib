/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef _MESHINTERFACEXHM_
#define _MESHINTERFACEXHM_

#include <memory>

#include "linearalgebra3d.h"
#include "mAOMD.h"
#include "mEdge.h"
#include "mFace.h"
#include "mRegion.h"
#include "mVertex.h"
#include "xMesh.h"
#include "xRegion.h"

namespace xfastmarching
{
/*
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
*/
#ifndef _MESHINTERFACEXREGION2_
// durty !!!!!!!!!!!!!! just to get filterVertex !!!!!!!!!!!!!!
#include "meshinterfacexRegion.h"
#endif

class meshinterfacexHM
{
  public:
   typedef AOMD::mVertex vertex;
   typedef AOMD::mEdge edge;
   typedef AOMD::mFace face;
   typedef AOMD::mRegion region;
   typedef AOMD::mEntity entity;
   typedef xfem::partmanxMesh_t partman_t;
   typedef xinterface::aomd::xAttachedDataManagerAOMD<bool>::c_iterKey_t side_iterator_t;
   /// iterator class used by this class
   class Iterator : public std::iterator<std::forward_iterator_tag, AOMD::mEntity *>
   //, ptrdiff_t, const AOMD::mEntity **, const AOMD::mEntity *&>
   {
     public:
      void findNext()
      {
         side_iterator_t end = HM->side[0].endKey();
         while (i != end)
         {
            auto &e = **i;
            if (!HM->background_to_iso.getData(e)) break;
            ++i;
         }
      }
      Iterator(const Iterator &rhs) : i(rhs.i), what(rhs.what), HM(rhs.HM) {}
      Iterator(const side_iterator_t &i_, int what_, const meshinterfacexHM *HM_) : i(i_), what(what_), HM(HM_)
      {
         if (!what) findNext();
      }
      reference operator*() { return const_cast<reference>(*i); }
      value_type operator*() const { return (*i); }
      Iterator &operator++()
      {
         assert(HM);
         ++i;
         if (what) return *this;
         findNext();
         return *this;
      }
      Iterator operator++(int)
      {
         Iterator tmp(*this);
         operator++();
         return tmp;
      }
      bool operator==(const Iterator &rhs) { return i == rhs.i; }
      bool operator!=(const Iterator &rhs) { return i != rhs.i; }

     private:
      side_iterator_t i;
      int what;
      const meshinterfacexHM *HM;
   };

   meshinterfacexHM(const xfem::xMesh &_iso_support_mesh, const xfem::xMesh &_background_mesh,
                    xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity *> &_background_to_iso,
                    xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity *> &_iso_to_background,
                    std::array<xinterface::aomd::xAttachedDataManagerAOMD<bool>, 4> &_side)
       : iso_support_mesh(_iso_support_mesh),
         background_mesh(_background_mesh),
         background_to_iso(_background_to_iso),
         iso_to_background(_iso_to_background),
         side(_side),
         f_reg(nullptr)
   {
   }
   Iterator begin(int what) const { return Iterator(side[what].beginKey(), what, this); }
   Iterator end(int what) const { return Iterator(side[what].endKey(), what, this); }
   typedef xfem::keyManagerSendOrReceive key_manager_sor_t;
   typedef xfem::keyManagerSendOrReceive::information_key_t information_key_t;
   typedef xfem::xFilteredRegion<partman_t::c_iter_object_t, filterVertex> filtered_region_t;
   typedef filtered_region_t::FilterIter FilterIter;
   std::unique_ptr<filtered_region_t> f_reg;
   const xfem::xMesh &iso_support_mesh;
   const xfem::xMesh &background_mesh;
   std::array<xinterface::aomd::xAttachedDataManagerAOMD<bool>, 4> &side;
   xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity *> &background_to_iso;
   xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity *> &iso_to_background;
};

inline std::size_t getTag(const meshinterfacexHM &mi, const std::string &tagname)
{
   return AOMD::AOMD_Util::Instance()->lookupMeshDataId(tagname.c_str());
};

inline std::size_t getNewTag(const meshinterfacexHM &mi, const std::string &tagname)
{
   return AOMD::AOMD_Util::Instance()->newMeshDataId(tagname.c_str());
}

inline void releaseTag(const meshinterfacexHM &mi, const std::string &tagname)
{
   int tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId(tagname.c_str());
   AOMD::AOMD_Util::Instance()->deleteMeshDataId(tag);
}

inline void releaseTag(const meshinterfacexHM &mi, const std::size_t &tag) { AOMD::AOMD_Util::Instance()->deleteMeshDataId(tag); }

#ifndef _MESHINTERFACEXREGION2_
class AttachableDataPointer : public AOMD::mAttachableData
{
  public:
   AttachableDataPointer(void *pnt_) : pnt(pnt_) {}
   ~AttachableDataPointer() override = default;
   void *pnt;
};
#endif

inline void *getAttachedDataPointer(const meshinterfacexHM &mi, AOMD::mEntity &en, size_t tag)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(en.getData(tag));
   if (dataatt)
      return dataatt->pnt;
   else
      return nullptr;
}

inline void deleteData(const meshinterfacexHM &mi, AOMD::mEntity &en, size_t tag) { en.deleteData(tag); }
inline void attachDataPointer(const meshinterfacexHM &mi, AOMD::mEntity &en, size_t tag, void *data)
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

inline void *getAttachedDataPointer(const meshinterfacexHM &mi, AOMD::mVertex &v, size_t tag)
{
   AttachableDataPointer *dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
      return dataatt->pnt;
   else
      return nullptr;
}

inline void attachDataPointer(const meshinterfacexHM &mi, AOMD::mVertex &v, size_t tag, void *data)
{
   AttachableDataPointer *dataatt;
   auto **pv = mi.background_to_iso.getData(static_cast<AOMD::mEntity &>(v));
   if (pv)
      dataatt = static_cast<AttachableDataPointer *>((*pv)->getData(tag));
   else
      dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
   {
      dataatt->pnt = data;
      return;
   }
   dataatt = new AttachableDataPointer(data);
   if (pv)
      (*pv)->attachData(tag, dataatt);
   else
      v.attachData(tag, dataatt);
}

inline void deleteData(const meshinterfacexHM &mi, AOMD::mVertex &v, size_t tag)
{
   auto **pv = mi.background_to_iso.getData(static_cast<AOMD::mEntity &>(v));
   if (pv)
      (*pv)->deleteData(tag);
   else
      v.deleteData(tag);
}

inline void *getAttachedDataPointer(const meshinterfacexHM &mi, AOMD::mFace &v, size_t tag)
{
   AttachableDataPointer *dataatt;
   auto **pv = mi.background_to_iso.getData(static_cast<AOMD::mEntity &>(v));
   if (pv)
      dataatt = static_cast<AttachableDataPointer *>((*pv)->getData(tag));
   else
      dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
      return dataatt->pnt;
   else
      return nullptr;
}

inline void attachDataPointer(const meshinterfacexHM &mi, AOMD::mFace &v, size_t tag, void *data)
{
   AttachableDataPointer *dataatt;
   auto **pv = mi.background_to_iso.getData(static_cast<AOMD::mEntity &>(v));
   if (pv)
      dataatt = static_cast<AttachableDataPointer *>((*pv)->getData(tag));
   else
      dataatt = static_cast<AttachableDataPointer *>(v.getData(tag));
   if (dataatt)
   {
      dataatt->pnt = data;
      return;
   }
   dataatt = new AttachableDataPointer(data);
   v.attachData(tag, dataatt);
}

inline void deleteData(const meshinterfacexHM &mi, AOMD::mFace &v, size_t tag) { v.deleteData(tag); }

template <class VERTEXINPUTITERATOR>
void getVerticesNeighbors(const meshinterfacexHM &mi, const AOMD::mVertex &v, VERTEXINPUTITERATOR it)
{
   size_t nbedges = v.size(1);
   for (size_t i = 0; i < nbedges; ++i)
   {
      AOMD::mEntity *e = v.get(1, i);
      if (mi.side[1].getData(*e))
      {
         for (size_t k = 0; k < 2; ++k)
         {
            auto *vike = e->get(0, k);
            const AOMD::mVertex *vik = static_cast<const AOMD::mVertex *>(vike);
            if (vik != &v)
            {
               auto **pv = mi.background_to_iso.getData(*vike);
               if (pv)
                  *it++ = static_cast<const AOMD::mVertex *>(*pv);
               else
                  *it++ = vik;
            }
         }
      }
   }
   auto **pv = mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v));
   if (!pv) pv = mi.iso_to_background.getData(static_cast<const AOMD::mEntity &>(v));
   if (pv)
   {
      auto &vi = **pv;
      nbedges = vi.size(1);
      for (size_t i = 0; i < nbedges; ++i)
      {
         AOMD::mEntity *e = vi.get(1, i);
         if (mi.side[1].getData(*e))
         {
            for (size_t k = 0; k < 2; ++k)
            {
               auto *vike = e->get(0, k);
               if (vike != *pv)
               {
                  // assert(!mi.background_to_iso.getData(*vike));
                  const AOMD::mVertex *vik = static_cast<const AOMD::mVertex *>(vike);
                  *it++ = vik;
               }
            }
         }
      }
   }
}

template <class EDGEINPUTITERATOR>
void getEdges(const meshinterfacexHM &mi, const AOMD::mVertex &v, EDGEINPUTITERATOR it)
{
   size_t n = v.size(1);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(1, i);
      if (mi.side[1].getData(*e))
      {
         *it = static_cast<AOMD::mEdge *>(e);
         ++it;
      }
   }
   auto **pv = mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v));
   if (!pv) pv = mi.iso_to_background.getData(static_cast<const AOMD::mEntity &>(v));
   if (pv)
   {
      auto &vi = **pv;
      n = vi.size(1);
      for (size_t i = 0; i < n; ++i)
      {
         auto *e = vi.get(1, i);
         if (mi.side[1].getData(*e))
         {
            *it = static_cast<AOMD::mEdge *>(e);
            ++it;
         }
      }
   }
}

template <class FACEINPUTITERATOR>
void getFaces(const meshinterfacexHM &mi, const AOMD::mVertex &v, FACEINPUTITERATOR it)
{
   size_t n = v.size(2);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(2, i);
      if (mi.side[2].getData(*e))
      {
         *it = static_cast<AOMD::mFace *>(e);
         ++it;
      }
   }
   auto **pv = mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v));
   if (!pv) pv = mi.iso_to_background.getData(static_cast<const AOMD::mEntity &>(v));
   if (pv)
   {
      auto &vi = **pv;
      n = vi.size(2);
      for (size_t i = 0; i < n; ++i)
      {
         auto *e = vi.get(2, i);
         if (mi.side[2].getData(*e))
         {
            *it = static_cast<AOMD::mFace *>(e);
            ++it;
         }
      }
   }
}

template <class REGIONINPUTITERATOR>
void getRegions(const meshinterfacexHM &mi, const AOMD::mVertex &v, REGIONINPUTITERATOR it)
{
   size_t n = v.size(3);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(3, i);
      if (mi.side[3].getData(*e))
      {
         *it = static_cast<AOMD::mRegion *>(e);
         ++it;
      }
   }
   auto **pv = mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v));
   if (!pv) pv = mi.iso_to_background.getData(static_cast<const AOMD::mEntity &>(v));
   if (pv)
   {
      auto &vi = **pv;
      n = vi.size(3);
      for (size_t i = 0; i < n; ++i)
      {
         auto *e = vi.get(3, i);
         if (mi.side[3].getData(*e))
         {
            *it = static_cast<AOMD::mRegion *>(e);
            ++it;
         }
      }
   }
}

template <class FACEINPUTITERATOR>
void getFaces(const meshinterfacexHM &mi, const AOMD::mEdge &v, FACEINPUTITERATOR it)
{
   size_t n = v.size(2);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(2, i);
      if (mi.side[2].getData(*e))
      {
         *it = static_cast<AOMD::mFace *>(e);
         ++it;
      }
   }
   auto **pv = mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v));
   if (!pv) pv = mi.iso_to_background.getData(static_cast<const AOMD::mEntity &>(v));
   if (pv)
   {
      auto &vi = **pv;
      n = vi.size(2);
      for (size_t i = 0; i < n; ++i)
      {
         auto *e = vi.get(2, i);
         if (mi.side[2].getData(*e))
         {
            *it = static_cast<AOMD::mFace *>(e);
            ++it;
         }
      }
   }
}

template <class REGIONINPUTITERATOR>
void getRegions(const meshinterfacexHM &mi, const AOMD::mEdge &v, REGIONINPUTITERATOR it)
{
   size_t n = v.size(3);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(3, i);
      if (mi.side[3].getData(*e))
      {
         *it = static_cast<AOMD::mRegion *>(e);
         ++it;
      }
   }
   auto **pv = mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v));
   if (!pv) pv = mi.iso_to_background.getData(static_cast<const AOMD::mEntity &>(v));
   if (pv)
   {
      auto &vi = **pv;
      size_t n = vi.size(3);
      for (size_t i = 0; i < n; ++i)
      {
         auto *e = vi.get(3, i);
         if (mi.side[3].getData(*e))
         {
            *it = static_cast<AOMD::mRegion *>(e);
            ++it;
         }
      }
   }
}

template <class VERTEXINPUTITERATOR>
void getVertices(const meshinterfacexHM &mi, const AOMD::mEdge &v, VERTEXINPUTITERATOR it)
{
   // assert(mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v)));
   size_t n = v.size(0);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(0, i);
      auto **pv = mi.background_to_iso.getData(*e);
      if (pv) e = *pv;
      if (mi.side[0].getData(*e))
      {
         *it = static_cast<AOMD::mVertex *>(e);
         ++it;
      }
   }
}

template <class VERTEXINPUTITERATOR>
void getVertices(const meshinterfacexHM &mi, const AOMD::mFace &v, VERTEXINPUTITERATOR it)
{
   // assert(mi.background_to_iso.getData(static_cast<AOMD::mEntity &>(v)));
   size_t n = v.size(0);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(0, i);
      auto **pv = mi.background_to_iso.getData(*e);
      if (pv) e = *pv;
      if (mi.side[0].getData(*e))
      {
         *it = static_cast<AOMD::mVertex *>(e);
         ++it;
      }
   }
}

template <class VERTEXINPUTITERATOR>
void getVertices(const meshinterfacexHM &mi, const AOMD::mRegion &v, VERTEXINPUTITERATOR it)
{
   // assert(mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v)));
   size_t n = v.size(0);
   for (size_t i = 0; i < n; ++i)
   {
      auto *e = v.get(0, i);
      auto **pv = mi.background_to_iso.getData(*e);
      if (pv) e = *pv;
      if (mi.side[0].getData(*e))
      {
         *it = static_cast<AOMD::mVertex *>(e);
         ++it;
      }
   }
}

inline double getX(const meshinterfacexHM &mi, const AOMD::mVertex &v) { return v.point()(0); }

inline double getY(const meshinterfacexHM &mi, const AOMD::mVertex &v) { return v.point()(1); }

inline double getZ(const meshinterfacexHM &mi, const AOMD::mVertex &v) { return v.point()(2); }

inline int getId(const meshinterfacexHM &mi, const AOMD::mVertex &v)
{
   auto **pv = mi.background_to_iso.getData(static_cast<const AOMD::mEntity &>(v));
   if (pv)
      return (*pv)->getId();
   else
      return v.getId();
}

template <class GEOMVECT>
void getCoord(const meshinterfacexHM &mi, const AOMD::mVertex &v, GEOMVECT &x)
{
   std::cout << "getcoord not implemented for the used type " << std::endl;
   assert(0);
}

template <class T>
void getCoord(const meshinterfacexHM &mi, const AOMD::mVertex &v, vector2d<T> &x)
{
   x = {v.point()(0), v.point()(1)};
}

template <class T>
void getCoord(const meshinterfacexHM &mi, const AOMD::mVertex &v, vector3d<T> &x)
{
   x = {v.point()(0), v.point()(1), v.point()(2)};
}

/*
inline int getNbVertex(const meshinterfacexHM &mi) { return mi.reg.size(0); }

inline const AOMD::mVertex *getRemoteVertex(const meshinterfacexHM &mi, const AOMD::mVertex &v, int remot_proc_id)
{
   auto po = mi.reg.getPartitionManager().getConstPartitionObject(static_cast<const AOMD::mEntity &>(v));
   return static_cast<const AOMD::mVertex *>(po.getRemoteObjectOn(remot_proc_id));
}

inline void setFilteredRegion(const meshinterfacexHM &mi)
{
   throw -458;
   return;
}

inline citerobjecttovertexiteratorconvertor<meshinterfacexHM::FilterIter> beginBnd(const meshinterfacexHM &mi)
{
   return citerobjecttovertexiteratorconvertor<meshinterfacexHM::FilterIter>(mi.f_reg->begin());
}
inline citerobjecttovertexiteratorconvertor<meshinterfacexHM::FilterIter> endBnd(const meshinterfacexHM &mi)
{
   return citerobjecttovertexiteratorconvertor<meshinterfacexHM::FilterIter>(mi.f_reg->end());
}
*/

}  // namespace xfastmarching
#endif
