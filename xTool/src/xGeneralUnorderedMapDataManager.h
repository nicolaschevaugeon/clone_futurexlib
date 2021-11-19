/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _xGeneralUnorderedMapDataManager_
#define _xGeneralUnorderedMapDataManager_
#include <iterator>
#include <unordered_map>

#include "xIteratorTools.h"

namespace xtool
{
/// \brief This class is a template implementation of an abstract DATAMANAGER interface (or concept) dedicated to key that
//! do not exist as an unique object n memory but as a instance that have a intrinsec unique representation.
///
///  xGeneralUnorderedMapDataManager is a simple implementation of the DATAMANAGER using an unordered map to store the relations.
///  xAttachedDataMnagerAOMD and xUnorderedMapDataManager are other implementation of the DATAMNANAGER.
/// Compare to xUnorderedMapDataManager the key used in the underlying associative container is not a pointer to it but an
/// instance of it. Some how it breaks the concept of attaching information to a concrete object given by its address but the
/// concept of attaching an information to a fictitious object remain. A  key define a unique fictitious object that will uniquely
/// be identified by the xGeneralUnorderedMapDataManager not mater that used concrete object are different.
//
template <typename KEYTYPE, typename DATATYPE, typename HASH_KEYTYPE = std::hash<KEYTYPE>,
          typename EQUAL_KEYTYPE = std::equal_to<KEYTYPE>>
class xGeneralUnorderedMapDataManager
{
  private:
   typedef std::unordered_map<KEYTYPE, DATATYPE *, HASH_KEYTYPE, EQUAL_KEYTYPE> data_container_t;
   data_container_t data_map;

  public:
   using value_type = DATATYPE;
   using key_type = KEYTYPE;

   /// Default constructor
   inline xGeneralUnorderedMapDataManager() = default;
   /// Copy constructor
   inline xGeneralUnorderedMapDataManager(const xGeneralUnorderedMapDataManager &in)
   {
      for (const auto &key_data : in.data_map)
      {
         data_map[key_data.first] = new DATATYPE(*(key_data.second));
      }
   }
   /// create a new data associated to key, and return a reference to it. If a data associated to key already exist in the
   /// xGeneralUnorderedMapDataManager instance, return a reference to the existing data. (Note: the semantic is the same as the
   /// operator[] of stl's associative container )
   inline DATATYPE &setData(const KEYTYPE &key)
   {
      DATATYPE *data = getData(key);
      if (!data)
      {
         data = new DATATYPE;
         data_map[key] = data;
      }
      return *data;
   }
   /*! return reference to object DATATYPE associated to e. If already exist return the existing one,
       otherwise construct a DATATYPE  by calling new DATATYPE(a) and associate it to e, before returning the ref
   */
   template <class... A>
   inline DATATYPE &emplaceSetData(const KEYTYPE &key, A &&... a)
   {
      DATATYPE *data = getData(key);
      if (!data)
      {
         data = new DATATYPE(std::forward<A>(a)...);
         data_map[key] = data;
      }
      return *data;
   }
   /// return a const pointer to the data associated to key
   inline const DATATYPE *getData(const KEYTYPE &key) const
   {
      const DATATYPE *data = nullptr;
      auto it = data_map.find(key);
      if (it != data_map.end()) data = it->second;
      return data;
   }
   /// return a pointer to the data associated to key
   inline DATATYPE *getData(const KEYTYPE &key)
   {
      DATATYPE *data = nullptr;
      auto it = data_map.find(key);
      if (it != data_map.end()) data = it->second;
      return data;
   }
   /// return a const ref to the data associated to key or throw if there is no data associated to key
   inline const DATATYPE &at(const KEYTYPE &key) const { return *data_map.at(const_cast<KEYTYPE *>(&key)); }
   /// return a ref to the data associated to key or throw if there is no data associated to key
   inline DATATYPE &at(const KEYTYPE &key) { return *data_map.at(&key); }
   /// return true if the key exist in the DataManager
   inline bool contains(const KEYTYPE &key) const
   {
      return data_map.count(&key);
      //! note : could use .contains in c++20.
   }
   /// delete data associated to key. if no data is associated to key, deleteData does nothing.
   /// all references to the data associated to key keeped by the user become invalide
   inline void deleteData(const KEYTYPE &key)
   {
      auto it = data_map.find(key);
      if (it != data_map.end())
      {
         delete it->second;
         data_map.erase(it);
      }
   }
   /// remove all delete all the data and the data/keys association (clear the data_map)
   inline void clear()
   {
      for (auto e : data_map)
      {
         delete e.second;
      }
      data_map.clear();
   }
   /// destructor :  clear the data map.
   inline ~xGeneralUnorderedMapDataManager() { clear(); }
   // iterate on data map keys is offered by use of special iterator
   /// Iterator on keys
   typedef std::iterator<std::forward_iterator_tag, typename data_container_t::key_type, std::ptrdiff_t,
                         const typename data_container_t::key_type *, const typename data_container_t::key_type &>
       std_iterator_t;

   class constIterKey : public std_iterator_t
   {
     public:
      constIterKey() = default;
      constIterKey(const constIterKey &rhs) : it(rhs.it) {}
      constIterKey &operator=(const constIterKey &rhs)
      {
         it = rhs.it;
         return *this;
      }
      typename std_iterator_t::reference operator*() const { return it->first; }
      constIterKey &operator++()
      {
         ++it;
         return *this;
      }
      constIterKey operator++(int)
      {
         constIterKey tmp(*this);
         ++it;
         return tmp;
      }
      bool operator==(const constIterKey &rhs) { return it == rhs.it; }
      bool operator!=(const constIterKey &rhs) { return it != rhs.it; }
      friend class xGeneralUnorderedMapDataManager;

     private:
      constIterKey(const typename data_container_t::const_iterator &it_) : it(it_) {}
      typename data_container_t::const_iterator it;
   };
   typedef constIterKey c_iterKey_t;

   inline c_iterKey_t beginKey() const
   {
      auto it = data_map.begin();
      return c_iterKey_t(it);
   }

   inline c_iterKey_t endKey() const
   {
      auto it = data_map.end();
      return c_iterKey_t(it);
   }
   inline xRange<c_iterKey_t> rangeKey() const { return xRange<c_iterKey_t>(beginKey(), endKey()); }
   inline bool empty() const { return data_map.empty(); }
   inline size_t size() const { return data_map.size(); }
};

}  // namespace xtool

#endif
