#ifndef _xAttachedDataManagerAOMD_
#define _xAttachedDataManagerAOMD_
#include <chrono>
#include <iostream>
#include <list>
#include <sstream>

#include "mAOMD.h"
#include "mAttachableDataContainer.h"
#include "mEntity.h"
#include "xIteratorTools.h"

namespace xinterface
{
namespace aomd
{
/// \brief This class is an implementation of abstract DATAMANAGER concept proposed in xUnorderedMapDataManager.
///
///  It treat the case of KEYTYPE being  a  mAttachableDataContainer type (i.e. all derived class mEntity, mFace,.... are
///  possible) DATATYPE remain a template parameter on the type of data to store.

template <typename KEYTYPE, typename DATATYPE>
class xAttachedDataManagerAOMDGeneric
{
  public:
   using value_type = DATATYPE;
   using key_type = KEYTYPE;
   typedef std::list<KEYTYPE *> key_container_t;
   typedef typename key_container_t::const_iterator c_iterKey_t;
   inline xAttachedDataManagerAOMDGeneric();
   inline xAttachedDataManagerAOMDGeneric(const xAttachedDataManagerAOMDGeneric &in);
   inline xAttachedDataManagerAOMDGeneric(xAttachedDataManagerAOMDGeneric &&in);
   inline ~xAttachedDataManagerAOMDGeneric();
   inline xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE> &operator=(xAttachedDataManagerAOMDGeneric in);

   inline const DATATYPE *getData(const KEYTYPE &e) const;
   inline DATATYPE *getData(const KEYTYPE &e);
   inline const DATATYPE &at(const KEYTYPE &key) const;
   inline DATATYPE &at(const KEYTYPE &key);
   /// return true if the key exist in the DataManager
   inline bool contains(const KEYTYPE &key) const;
   /*! return reference to object DATATYPE associated to e. If already exist return the existing one,
       otherwise default construct a DATATYPE and associate it to e, before returning the ref
   */
   inline DATATYPE &setData(const KEYTYPE &e);
   /*! return reference to object DATATYPE associated to e. If already exist return the existing one,
       otherwise construct a DATATYPE  by calling new DATATYPE(a) and associate it to e, before returning the ref
   */
   template <class... A>
   inline DATATYPE &emplaceSetData(const KEYTYPE &e, A &&... a);
   /// try to delete data attached to entity e.
   /*! return 1 if data deleted, 0 otherwise !*/
   inline size_t deleteData(const KEYTYPE &e);
   inline void clear();
   inline void swap(xAttachedDataManagerAOMDGeneric &in);

   inline c_iterKey_t beginKey() const;
   inline c_iterKey_t endKey() const;
   /// return the range of keys known by the container
   inline xtool::xRange<c_iterKey_t> rangeKey() const;

   inline bool empty() const { return garbage.empty(); }
   inline size_t size() const { return garbage.size(); }

  protected:
   typedef typename key_container_t::iterator iterKey_t;
   size_t tag;
   bool zombi;
   std::list<KEYTYPE *> garbage;
   class xAttachablePnt : public AOMD::mAttachableData
   {
     public:
      DATATYPE *pnt = nullptr;
      iterKey_t it;
   };
   void setTag();
};

//! original API concerved by use of this derivation.
// Hoppefuly it will disapear one day ...
// So do not use it if you can
template <class DATATYPE>
class xAttachedDataManagerAOMD : public xAttachedDataManagerAOMDGeneric<AOMD::mEntity, DATATYPE>
{
   // please please don't ever think about adding things here
};

// =====================================================
// implementation
// =====================================================
template <typename KEYTYPE, typename DATATYPE>
inline xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::xAttachedDataManagerAOMDGeneric() : tag(0), zombi(false)
{
   static_assert(std::is_base_of<AOMD::mEntity, KEYTYPE>::value, "KEYTYPE must be derived from AOMD::mEntity");
   setTag();
}

template <typename KEYTYPE, typename DATATYPE>
inline xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::xAttachedDataManagerAOMDGeneric(xAttachedDataManagerAOMDGeneric &&in)
    : tag(in.tag), zombi(false), garbage(std::move(in.garbage))
{
   in.zombi = true;
   // std::cout << "Move Constructor, tag : "<< tag << std::endl;
}

template <typename KEYTYPE, typename DATATYPE>
inline void xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::swap(xAttachedDataManagerAOMDGeneric &in)
{
   std::swap(tag, in.tag);
   std::swap(zombi, in.zombi);
   std::swap(garbage, in.garbage);
}

template <typename KEYTYPE, typename DATATYPE>
inline xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE> &xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::operator=(
    xAttachedDataManagerAOMDGeneric in)
{
   swap(in);
   return *this;
}

template <typename KEYTYPE, typename DATATYPE>
inline xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::xAttachedDataManagerAOMDGeneric(
    const xAttachedDataManagerAOMDGeneric &in)
    : tag(0), zombi(false)
{
   setTag();
   for (auto e : in.garbage)
   {
      auto data = in.getData((*e));
      if (data)
      {
         setData(*e) = *data;
      }
   }
   // std::cout << "Copy Constructor, tag :"<< tag << std::endl;
}

template <typename KEYTYPE, typename DATATYPE>
inline const DATATYPE *xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::getData(const KEYTYPE &e) const
{
   auto attachable_pnt = static_cast<xAttachablePnt *>(e.getData(tag));
   if (attachable_pnt) return attachable_pnt->pnt;
   return nullptr;
}

template <typename KEYTYPE, typename DATATYPE>
inline DATATYPE *xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::getData(const KEYTYPE &e)
{
   auto attachable_pnt = static_cast<xAttachablePnt *>(e.getData(tag));
   if (attachable_pnt) return attachable_pnt->pnt;
   return nullptr;
}

template <typename KEYTYPE, typename DATATYPE>
inline bool xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::contains(const KEYTYPE &key) const
{
   return key.getData(tag);
}

/// return a const ref to the data associated to key or throw if there is no data associated to key
template <typename KEYTYPE, typename DATATYPE>
inline const DATATYPE &xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::at(const KEYTYPE &key) const
{
   const DATATYPE *data = getData(key);
   if (!data) throw;
   return *data;
}

/// return a ref to the data associated to key or throw if there is no data associated to key
template <typename KEYTYPE, typename DATATYPE>
inline DATATYPE &xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::at(const KEYTYPE &key)
{
   DATATYPE *data = getData(key);
   if (!data) throw;
   return *data;
}

template <typename KEYTYPE, typename DATATYPE>
inline DATATYPE &xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::setData(const KEYTYPE &e)
{
   // std::cout << "xAttachedDataManagerAOMDGeneric setData TAG " << tag << " e " << &e << std::endl;
   DATATYPE *data = getData(e);
   if (!data)
   {
      KEYTYPE &e_nc = const_cast<KEYTYPE &>(e);
      data = new DATATYPE;
      xAttachablePnt *adata = new xAttachablePnt;
      adata->pnt = data;
      e_nc.attachData(tag, adata);
      // adata->it = garbage.insert(garbage.begin(), &e);
      garbage.push_back(&e_nc);
      adata->it = --(garbage.end());
   }
   return *data;
}

template <typename KEYTYPE, typename DATATYPE>
template <class... A>
inline DATATYPE &xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::emplaceSetData(const KEYTYPE &e, A &&... a)
{
   DATATYPE *data = getData(e);
   if (!data)
   {
      KEYTYPE &e_nc = const_cast<KEYTYPE &>(e);
      data = new DATATYPE(std::forward<A>(a)...);
      xAttachablePnt *adata = new xAttachablePnt;
      adata->pnt = data;
      e_nc.attachData(tag, adata);
      // adata->it = garbage.insert(garbage.begin(), &e);
      garbage.push_back(&e_nc);
      adata->it = --(garbage.end());
   }
   return *data;
}

template <typename KEYTYPE, typename DATATYPE>
inline size_t xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::deleteData(const KEYTYPE &e)
{
   auto attachable_pnt = static_cast<xAttachablePnt *>(e.getData(tag));
   if (attachable_pnt)
   {
      //	std::cout << "IN attachable Data deleteData on "<< &e <<  " garbage size "<< garbage.size() << std::endl;
      delete (attachable_pnt->pnt);
      garbage.erase(attachable_pnt->it);
      const_cast<KEYTYPE &>(e).deleteData(tag);
      return 1;
      //	std::cout << "Done attachable Data deleteData garbage size " << garbage.size() << std::endl;
   }
   return 0;
}

template <typename KEYTYPE, typename DATATYPE>
inline void xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::clear()
{
   for (auto it = garbage.begin(), ite = garbage.end(); it != ite;)
   {
      KEYTYPE *e = *it;
      ++it;
      deleteData(*e);
   }
}

template <typename KEYTYPE, typename DATATYPE>
inline xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::~xAttachedDataManagerAOMDGeneric()
{
   if (!zombi)
   {
      //  std::cout << "destructor, tag :" << tag << std::endl;
      clear();
      AOMD::AOMD_Util::Instance()->deleteMeshDataId(tag);
   }
}

template <typename KEYTYPE, typename DATATYPE>
inline auto xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::beginKey() const -> c_iterKey_t
{
   if (zombi)
      throw -1;
   else
      return garbage.cbegin();
}

template <typename KEYTYPE, typename DATATYPE>
inline auto xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::endKey() const -> c_iterKey_t
{
   if (zombi)
      throw -1;
   else
      return garbage.cend();
}

template <typename KEYTYPE, typename DATATYPE>
inline auto xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::rangeKey() const -> xtool::xRange<c_iterKey_t>
{
   return xtool::xRange<c_iterKey_t>(beginKey(), endKey());
}

template <typename KEYTYPE, typename DATATYPE>
inline void xAttachedDataManagerAOMDGeneric<KEYTYPE, DATATYPE>::setTag()
{
   std::stringstream tagname;
   std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
   std::chrono::high_resolution_clock::duration dt = t.time_since_epoch();
   tagname << "xAttachedDataManagerAOMDGeneric" << this << dt.count();
   tag = AOMD::AOMD_Util::Instance()->newMeshDataId(tagname.str().c_str());
}

}  // namespace aomd
}  // namespace xinterface

#endif
