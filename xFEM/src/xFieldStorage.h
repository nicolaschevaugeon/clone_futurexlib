/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XFIELDSTORAGE_H
#define _XFIELDSTORAGE_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

#include "xApproxFctPtr.h"
#include "xSpace.h"
#include "xValueManager.h"

namespace xfem
{
/// A storage policy for the pointer to the values in double manager.
/*!
  This is the base class that define the interface.
  to class are derived from it so far :
  xFieldStorageFull and xFieldStorageElement. Note that the strorage is not necessarly thread safe ....
  !*/
template <typename VT>
class xFieldStorage
{
  public:
   typedef spacePtr_t spacePtr;
   typedef approxFctPtr_t shapeFctPtr;
   typedef femFcts_t femFcts;
   typedef std::vector<xValKey> femKeys;

   xFieldStorage(xValueManagerDist<VT>& d, const std::vector<spacePtr>& s) : double_manager(d), spaces(s){};
   virtual ~xFieldStorage() = default;
   ;
   virtual void clear() = 0;
   virtual int sizeFcts(AOMD::mEntity* e) const = 0;
   virtual std::vector<shapeFctPtr>::iterator beginFcts(AOMD::mEntity* e) const = 0;
   virtual std::vector<shapeFctPtr>::iterator endFcts(AOMD::mEntity* e) const = 0;
   virtual typename std::vector<xValue<VT>*>::iterator beginValues(AOMD::mEntity* e) const = 0;
   virtual typename std::vector<xValue<VT>*>::iterator endValues(AOMD::mEntity* e) const = 0;

  protected:
   typedef typename std::vector<xValue<VT>*> values_t;
   typedef femFcts functions_t;
   mutable functions_t::iterator current_beginFcts;
   mutable functions_t::iterator current_endFcts;
   mutable typename values_t::iterator current_beginValues;
   mutable typename values_t::iterator current_endValues;
   xValueManagerDist<VT>& double_manager;
   const std::vector<spacePtr>& spaces;
};

/// A storage policy for the pointer to the values in double manager.
/*!
  A bigger storage for higher perf if value are accessed using the field multiple time
  (for example,  explicit problem where you compute n time the right hand side)
  Until it is cleared calling clear,  this will store all the functions and values associated to each element that you have asked
  for. It can become huge in memory ...
  !*/
template <typename VT>
class xFieldStorageFull : public xFieldStorage<VT>
{
  public:
   typedef spacePtr_t spacePtr;
   typedef approxFctPtr_t shapeFctPtr;
   xFieldStorageFull(xValueManagerDist<VT>& d, const std::vector<spacePtr>& s) : xFieldStorage<VT>(d, s), e_current(nullptr) {}
   void clear() override;
   int sizeFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::functions_t::iterator beginFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::functions_t::iterator endFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::values_t::iterator beginValues(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::values_t::iterator endValues(AOMD::mEntity* e) const override;

  protected:
   mutable AOMD::mEntity* e_current;
   typedef std::pair<typename xFieldStorage<VT>::values_t, typename xFieldStorage<VT>::functions_t> stored_data_t;
   typedef std::unordered_map<AOMD::mEntity*, stored_data_t> hash_map_t;
   typedef typename hash_map_t::const_iterator hash_map_const_iterator_t;
   typedef typename hash_map_t::iterator hash_map_iterator_t;
   mutable hash_map_t data_map;
   void getApproximation(AOMD::mEntity* e) const;
   hash_map_iterator_t storeApproximation(AOMD::mEntity* e) const;
};

/// A storage policy for the pointer to the values in double manager.
/*!
  Only values and approxfunction for the "current" element are stored.
  This is the default storage used in xField.
  !*/
template <typename VT>
class xFieldStorageElement : public xFieldStorage<VT>
{
  public:
   typedef spacePtr_t spacePtr;
   typedef approxFctPtr_t shapeFctPtr;
   xFieldStorageElement(xValueManagerDist<VT>& d, const std::vector<spacePtr>& s) : xFieldStorage<VT>(d, s), e_current(nullptr) {}
   void clear() override;
   int sizeFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::functions_t::iterator beginFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::functions_t::iterator endFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::values_t::iterator beginValues(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::values_t::iterator endValues(AOMD::mEntity* e) const override;

  protected:
   mutable AOMD::mEntity* e_current;
   mutable typename xFieldStorage<VT>::functions_t functions_current;
   mutable typename xFieldStorage<VT>::values_t values_current;
   void getApproximation(AOMD::mEntity* e) const;
};
/// A storage policy for the pointer to the values in double manager.
/*!
  Same kind as xFieldStorageFull except that race condition are not present as soon as all element are passed in review.
  During transition from empty to full it is not thread safe as realocation are possible during filling step.
  Until it is cleared calling clear,  this will store all the functions and values associated to each element that you have asked
  for. It can become huge in memory ...
  !*/
template <typename VT>
class xFieldStorageFullAttached : public xFieldStorage<VT>
{
  public:
   typedef spacePtr_t spacePtr;
   typedef approxFctPtr_t shapeFctPtr;
   typedef std::pair<size_t, size_t> data_t;

   xFieldStorageFullAttached(xValueManagerDist<VT>& d, const std::vector<spacePtr>& s) : xFieldStorage<VT>(d, s){};
   void clear() override;
   int sizeFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::functions_t::iterator beginFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::functions_t::iterator endFcts(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::values_t::iterator beginValues(AOMD::mEntity* e) const override;
   typename xFieldStorage<VT>::values_t::iterator endValues(AOMD::mEntity* e) const override;

  protected:
   mutable typename xFieldStorage<VT>::values_t values;
   mutable typename xFieldStorage<VT>::functions_t functions;
   mutable xinterface::aomd::xAttachedDataManagerAOMD<data_t> idxs;

   const data_t& getApproximation(AOMD::mEntity* e) const;
   const data_t* storeApproximation(AOMD::mEntity* e) const;
};

//--------------------------------------------------------------------------
template <typename VT>
int xFieldStorageFull<VT>::sizeFcts(AOMD::mEntity* e) const
{
   getApproximation(e);
   return xFieldStorage<VT>::current_endFcts - xFieldStorage<VT>::current_beginFcts;
};

template <typename VT>
typename xFieldStorage<VT>::functions_t::iterator xFieldStorageFull<VT>::beginFcts(AOMD::mEntity* e) const
{
   getApproximation(e);
   return xFieldStorage<VT>::current_beginFcts;
};

template <typename VT>
typename xFieldStorage<VT>::functions_t::iterator xFieldStorageFull<VT>::endFcts(AOMD::mEntity* e) const
{
   getApproximation(e);
   return xFieldStorage<VT>::current_endFcts;
};

template <typename VT>
typename xFieldStorage<VT>::values_t::iterator xFieldStorageFull<VT>::beginValues(AOMD::mEntity* e) const
{
   getApproximation(e);
   return xFieldStorage<VT>::current_beginValues;
};

template <typename VT>
typename xFieldStorage<VT>::values_t::iterator xFieldStorageFull<VT>::endValues(AOMD::mEntity* e) const
{
   getApproximation(e);
   return xFieldStorage<VT>::current_endValues;
};

template <typename VT>
void xFieldStorageFull<VT>::getApproximation(AOMD::mEntity* e) const
{
   if (e != e_current)
   {
      e_current = e;
      hash_map_iterator_t it = data_map.find(e_current);
      if (it == data_map.end())
      {
         it = storeApproximation(e_current);
      }

      xFieldStorage<VT>::current_beginValues = it->second.first.begin();
      xFieldStorage<VT>::current_endValues = it->second.first.end();
      xFieldStorage<VT>::current_beginFcts = it->second.second.begin();
      xFieldStorage<VT>::current_endFcts = it->second.second.end();
   }
}

template <typename VT>
typename xFieldStorageFull<VT>::hash_map_iterator_t xFieldStorageFull<VT>::storeApproximation(AOMD::mEntity* e) const
{
   // const bool debug = xdebug_flag;
   stored_data_t tmp;
   hash_map_iterator_t it = (data_map.insert(std::make_pair(e, tmp))).first;
   stored_data_t& data_e = it->second;
   typename xFieldStorage<VT>::values_t& values = data_e.first;
   typename xFieldStorage<VT>::functions_t& functions = data_e.second;
   values.clear();
   functions.clear();
#ifndef OLDVERSION_XFINITEELEMENT
   typename xFieldStorage<VT>::femKeys keys;
   for_each(xFieldStorage<VT>::spaces.begin(), xFieldStorage<VT>::spaces.end(),
            [&e, &keys, &functions](spacePtr f) { f->getKeysAndFcts(e, &keys, &functions); });
   values.reserve(functions.size());
   xFieldStorage<VT>::double_manager.getValPtr(keys.begin(), keys.end(), values);
#endif
#ifdef OLDVERSION_XFINITEELEMENT
   xFiniteElement FEM;
   FEM.setKeysAndFcts(e, spaces.begin(), spaces.end(), "trial");
   functions.assign(FEM.beginFct("trial"), FEM.endFct("trial"));
   xFieldStorage<VT>::double_manager.getValPtr(FEM.beginKey("trial"), FEM.endKey("trial"), values);
#endif
   return it;
}

template <typename VT>
void xFieldStorageFull<VT>::clear()
{
   e_current = nullptr;
   data_map.clear();
}

template <typename VT>
void xFieldStorageElement<VT>::clear()
{
   e_current = nullptr;
   values_current.clear();
   functions_current.clear();
}

template <typename VT>
void xFieldStorageElement<VT>::getApproximation(AOMD::mEntity* e) const
{
   if (e != e_current)
   {
      e_current = e;
      values_current.clear();
      functions_current.clear();
#ifndef OLDVERSION_XFINITEELEMENT
      typename xFieldStorage<VT>::femKeys keys;
      for_each(xFieldStorage<VT>::spaces.begin(), xFieldStorage<VT>::spaces.end(),
               [&e, &keys, this](spacePtr f) { f->getKeysAndFcts(e, &keys, &functions_current); });
      values_current.reserve(keys.size());
      xFieldStorage<VT>::double_manager.getValPtr(keys.begin(), keys.end(), values_current);
#endif
#ifdef OLDVERSION_XFINITEELEMENT
      xFiniteElement FEM;
      FEM.setKeysAndFcts(e, spaces.begin(), spaces.end(), "trial");
      functions_current.assign(FEM.beginFct("trial"), FEM.endFct("trial"));
      values_current.reserve(functions_current.size());
      xFieldStorage<VT>::double_manager.getValPtr(FEM.beginKey("trial"), FEM.endKey("trial"), values_current);
#endif
   }
}

template <typename VT>
int xFieldStorageElement<VT>::sizeFcts(AOMD::mEntity* e) const
{
   getApproximation(e);
   return functions_current.size();
};

template <typename VT>
typename xFieldStorage<VT>::functions_t::iterator xFieldStorageElement<VT>::beginFcts(AOMD::mEntity* e) const
{
   getApproximation(e);
   return functions_current.begin();
};

template <typename VT>
typename xFieldStorage<VT>::functions_t::iterator xFieldStorageElement<VT>::endFcts(AOMD::mEntity* e) const
{
   getApproximation(e);
   return functions_current.end();
};

template <typename VT>
typename xFieldStorage<VT>::values_t::iterator xFieldStorageElement<VT>::beginValues(AOMD::mEntity* e) const
{
   getApproximation(e);
   return values_current.begin();
};

template <typename VT>
typename xFieldStorage<VT>::values_t::iterator xFieldStorageElement<VT>::endValues(AOMD::mEntity* e) const
{
   getApproximation(e);
   return values_current.end();
};

//----------------------------------------
template <typename VT>
int xFieldStorageFullAttached<VT>::sizeFcts(AOMD::mEntity* e) const
{
   const data_t& r = getApproximation(e);
   return static_cast<int>(r.second - r.first);
}

template <typename VT>
typename xFieldStorage<VT>::functions_t::iterator xFieldStorageFullAttached<VT>::beginFcts(AOMD::mEntity* e) const
{
   return std::next(functions.begin(), getApproximation(e).first);
}

template <typename VT>
typename xFieldStorage<VT>::functions_t::iterator xFieldStorageFullAttached<VT>::endFcts(AOMD::mEntity* e) const
{
   return std::next(functions.begin(), getApproximation(e).second);
}

template <typename VT>
typename xFieldStorage<VT>::values_t::iterator xFieldStorageFullAttached<VT>::beginValues(AOMD::mEntity* e) const
{
   return std::next(values.begin(), getApproximation(e).first);
}

template <typename VT>
typename xFieldStorage<VT>::values_t::iterator xFieldStorageFullAttached<VT>::endValues(AOMD::mEntity* e) const
{
   return std::next(values.begin(), getApproximation(e).second);
}

template <typename VT>
const typename xFieldStorageFullAttached<VT>::data_t& xFieldStorageFullAttached<VT>::getApproximation(AOMD::mEntity* e) const
{
   const data_t* p = idxs.getData(*e);
   if (!p)
   {
      p = storeApproximation(e);
   }
   return *p;
}

template <typename VT>
const typename xFieldStorageFullAttached<VT>::data_t* xFieldStorageFullAttached<VT>::storeApproximation(AOMD::mEntity* e) const
{
   // const bool debug = xdebug_flag;
   data_t tmp;
   tmp.first = values.size();
   assert(tmp.first == functions.size());
   typename xFieldStorage<VT>::femKeys keys;
   auto& ref_functions = functions;
   for_each(xFieldStorage<VT>::spaces.begin(), xFieldStorage<VT>::spaces.end(),
            [&e, &keys, &ref_functions](spacePtr f) { f->getKeysAndFcts(e, &keys, &ref_functions); });
   values.reserve(functions.size());
   xFieldStorage<VT>::double_manager.getValPtr(keys.begin(), keys.end(), values);
   tmp.second = values.size();
   idxs.setData(*e) = tmp;
   return idxs.getData(*e);
}

template <typename VT>
void xFieldStorageFullAttached<VT>::clear()
{
   values.clear();
   functions.clear();
   idxs.clear();
}

//--------------------------------------------------------------------------

}  // namespace xfem
#endif  // _XFIELDSTORAGE_H
