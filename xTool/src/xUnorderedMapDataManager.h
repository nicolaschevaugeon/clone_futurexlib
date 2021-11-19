/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _xUnorderedMapDataManager_
#define _xUnorderedMapDataManager_
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <unordered_map>

#include "xIteratorTools.h"

namespace xtool
{
/// \brief This class is an implementation of an abstract DATAMANAGER interface (or concept).
///
///  DATAMANAGER is a class designed to store data (DATATYPE) relative to a key (KEYTYPE)
///  each DATAMANAGER implementation requires a default constructor, a copy constructor and the following member function  \n
///  DATATYPE & setData (const KEYTYPE &e) \n
///  DATATYPE &       emplaceSetData (const KEYTYPE &key,A &&...a )

///  const DATATYPE * getData( const KEYTYPE &e) const \n
///  DATATYPE * getData( const KEYTYPE &e) \ntemplate < class ...A >
///  const DATATYPE & at(const KEYTYPE &key) const;
///  DATATYPE & at(const KEYTYPE &key);
///  void   deleteData( const KEYTYPE &e ) \n
///  void clear() \n
///
///  One unique object of type DATATYPE exist and is owned by the DATAMANAGER for each key inserted to the DATAMANAGER using
///  setData \n DATATYPE & setData ( const KEYTYPE &e ) create new object in the DATAMANAGER, associated to key e., and returns a
///  reference to the new object \n deleteData( const KEYTYPE &e ) destroy the object associated to key e \n DATATYPE *
///  getData(KEYTYPE &e ) return a pointer to the object associated to key e if it exist, other wise it return a nullptr. \n
///  DATATYPE * getData(const KEYTYPE &e ) const return a pointer to the object associated to key e that is not modifyable. \n
///  void clear() empty the DATAMNAGER of all the key/objetc relation.\n
///  Note on the copy constructor : the DATAMANAGER return by the copy constructor, must have no relation left with the
///  DATAMANAGER it is copyied from it is supposed to be a deep copy ! \n
///
///  xUnorderedMapDataManager is a simple implementation of the DATAMANAGER using an unordered map to store the relations. It can
///  be used as an example to build other DATAMANAGER with the same concept. xAttachedDataMnagerAOMD is an other implementation of
///  the DATAMNANAGER that uses  mAttachedDataContainer object as keys and store the key/object relation (see
///  xinterface::aomd::xAttachedDataManagerAOMD.h) \n Note that this implementation store pointers to the keys. 2 keys are
///  considered equal if they have the same address. It means that a copy of a key is NOT the same key from the
///  xUnorderedMapDataManager point of view.
template <typename KEYTYPE, typename DATATYPE>
class xUnorderedMapDataManager
{
  private:
   typedef std::unordered_map<KEYTYPE *, DATATYPE *> data_container_t;
   data_container_t data_map;

  public:
   using value_type = DATATYPE;
   using key_type = KEYTYPE;

   /// Default constructor
   inline xUnorderedMapDataManager() = default;
   /// Copy constructor
   inline xUnorderedMapDataManager(const xUnorderedMapDataManager &in)
   {
      for (const auto &key_data : in.data_map)
      {
         data_map[key_data.first] = new DATATYPE(*(key_data.second));
      }
   }
   /// create a new data associated to key, and return a reference to it. If a data associated to key already exist in the
   /// xUnorderedMapDataManager instance, return a reference to the existing data. (Note: the semantic is the same as the
   /// operator[] of stl's associative container )
   inline DATATYPE &setData(const KEYTYPE &key)
   {
      DATATYPE *data = getData(key);
      if (!data)
      {
         data = new DATATYPE;
         data_map[const_cast<KEYTYPE *>(&key)] = data;
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
         data_map[const_cast<KEYTYPE *>(&key)] = data;
      }
      return *data;
   }

   /// return a const pointer to the data associated to key or nullptr if there is no data associated to key
   inline const DATATYPE *getData(const KEYTYPE &key) const
   {
      const DATATYPE *data = nullptr;
      //  it's seems ok for me to just pass &key to find() but g++ 4.8.4 does not agree ...
      auto it = data_map.find(const_cast<KEYTYPE *>(&key));
      if (it != data_map.end()) data = it->second;
      return data;
   }
   /// return a pointer to the data associated to key or nullptr if there is no data associated to key
   inline DATATYPE *getData(const KEYTYPE &key)
   {
      DATATYPE *data = nullptr;
      auto it = data_map.find(const_cast<KEYTYPE *>(&key));
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
      auto it = data_map.find(const_cast<KEYTYPE *>(&key));
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
   inline ~xUnorderedMapDataManager() { clear(); }
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
      friend class xUnorderedMapDataManager;

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
   /// return the range of keys known by the container
   inline xRange<c_iterKey_t> rangeKey() const { return xRange<c_iterKey_t>(beginKey(), endKey()); }
   inline bool empty() const { return data_map.empty(); }
   inline size_t size() const { return data_map.size(); }
};

}  // namespace xtool

#endif
