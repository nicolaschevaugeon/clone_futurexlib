/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef _xValKeyDataManager_H
#define _xValKeyDataManager_H
#include <iterator>
#include <unordered_map>
#include "xValKey.h"

namespace xfem
{

/// \brief This class is an implementation of an abstract DATAMANAGER interface (or concept) dedicated to  xValKey specific key
///
///  xValKeyDataManager is a simple implementation of the DATAMANAGER using an unordered map to store the relations.
///  xAttachedDataMnagerAOMD and xUnorderedMapDataManager are other implementation of the DATAMNANAGER.
/// Compare to xUnorderedMapDataManager the key used in the underlying associative container is a xValKey not a pointer to it.
/// Some how it breaks the concept of attaching information to a concrete object given by its address but the concept of
/// attaching an information to a fictitious object remain. A  xValKey define a unique fictitious object that will uniquely be
/// identified by the xValKeyDataManager not mater that used concrete object are different.
//
template < class DATATYPE >
class xValKeyDataManager
{
    private:
        typedef std::unordered_map < xValKey,  DATATYPE *,xHashValKey, xEqualValKey  > data_container_t;
        data_container_t data_map;
    public:

        /// Default constructor
        inline xValKeyDataManager() = default;
        /// Copy constructor
        inline xValKeyDataManager(const xValKeyDataManager & in)
        {
            for (const auto &key_data : in.data_map )
            {
                data_map[key_data.first] = new DATATYPE(*( key_data.second ));
            }
        }
        /// create a new data associated to key, and return a reference to it. If a data associated to key already exist in the xValKeyDataManager instance,
        /// return a reference to the existing data. (Note: the semantic is the same as the operator[] of stl's associative container )
        inline DATATYPE & setData (xValKey &key)
        {
            DATATYPE * data = getData(key);
            if (!data)
            {
                data = new DATATYPE;
                data_map[ key] = data;
            }
            return *data;
        }
        /// return a const pointer to the data associated to key
        inline const DATATYPE * getData( const xValKey &key) const
        {
            const DATATYPE *data = nullptr;
            auto it = data_map.find(key);
            if (it != data_map.end()) data = it->second;
            return data;
        }
        /// return a pointer to the data associated to key
        inline DATATYPE * getData(xValKey &key)
        {
            DATATYPE *data = nullptr;
            auto it = data_map.find(key);
            if (it != data_map.end()) data = it->second;
            return data;
        }
        /// delete data associated to key. if no data is associated to key, deleteData does nothing.
        /// all references to the data associated to key keeped by the user become invalide
        inline void   deleteData( xValKey &key )
        {
            auto it = data_map.find( key);
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
        inline ~xValKeyDataManager()
        {
            clear();
        }
        // iterate on data map keys is offered by use of special iterator
        /// Iterator on keys
        typedef std::iterator < std::forward_iterator_tag, typename data_container_t::key_type, std::ptrdiff_t, const typename data_container_t::key_type *, const typename data_container_t::key_type & >  std_iterator_t;

        class constIterKey : public std_iterator_t
        {
            public:
                constIterKey() = default;
                constIterKey(const constIterKey & rhs) : it(rhs.it){}
                constIterKey & operator= (const constIterKey & rhs) { it = rhs.it; return *this; }
                typename std_iterator_t::reference operator*() const { return it->first; }
                constIterKey & operator++() { ++it; return *this; }
                constIterKey operator++(int) { constIterKey tmp(*this); ++it; return tmp; }
                bool operator==(const constIterKey & rhs) { return it == rhs.it; }
                bool operator!=(const constIterKey & rhs) { return it != rhs.it; }
                friend class xValKeyDataManager;
            private:
                constIterKey(const typename data_container_t::const_iterator & it_) : it(it_) {}
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

};

} // end namespace

#endif
