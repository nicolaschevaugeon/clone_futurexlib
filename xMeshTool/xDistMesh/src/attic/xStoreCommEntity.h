/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef XSTORECOMMENTITY_H
#define XSTORECOMMENTITY_H

#include "xComm.h"

namespace distmesh
{

/// Template class model of storing container for xCommEntity informations corresponding to entity of any type
//! ET is the template parameter representing the type of entity we are working with
//! It is user responsibility to maintain a coherent use of process id between method of the class
template < typename ET >
class xStoreCommEntity
{
    public:
        /// Default constructor to initialize storage properly
        xStoreCommEntity(void);
        /// Copy constructor to initialize storage properly from an other instance
        xStoreCommEntity(const xStoreCommEntity& other);
        /// Destructor to clean correctly storage
        ~xStoreCommEntity(void);
        /// Method to add in the storage information relate to a duplicate entity in a remote process
        //! \param [in] ladress local address of an entity duplicated in a other process
        //! \param [in] rproc remote pid
        //! \param [in] radress remote address of the duplicate entity
        void set(ET *ladress, size_t rproc, ET *radress);
        /// Method to get xCommEntity information related to entity e
        //!  \return xCommEntity relate to e
        //!  \return null if no xCommEntity is associated to e
        xCommEntity * get(ET* e);
        /// Remove xCommEntity information related to e in the container
        void clear(ET* e);
        /// Remove  xCommConector information for rproc from xCommEntity related to e in the container
        //! If rproc was the last xCommConector in xCommEntity container then xCommEntity is removed
        //! from storage
        void clear(ET* e, size_t rproc);
        /// Assign a proc_id information to instance
        //! \param [in] proc_id_ id of the current process this instance belong to.
        void assignProcId(size_t proc_id_);
        /// Tells if e entity is the owner among all proc of the information associated to it
        //! This is used to create communication pattern where only owner are sending to all
        //! related proc. This is removing intermediate communication in between non owner proc
        //! \param [in] e entity address which is supposed to have information to communicate
        //! \return true if e is the owner (This include entity with no associated xCommEntity)
        //! \return false if e is not the owner and have an associated xCommEntity
        bool isOwner(ET *e);



        /*
         * TO BE DISCUSS
         * replaced by Exchanger and information familly..
         *
        // give location of array of ET addresses related to 
        // remote proc id rid and entity type what.
        // In those array address are sorted according to local or remote address value
        ET** getLocalAdressLocalSort(size_t rid, int what);
        ET** getLocalAdressRemotSort(size_t rid, int what);
        ET** getRemotAdressLocalSort(size_t rid, int what);
        ET** getRemotAdressRemotSort(size_t rid, int what);
        // same with owner stuff
        ET** getLocalAdressLocalSortOwnerOnly(size_t rid, int what);
        ET** getLocalAdressRemotSortOwnerOnly(size_t rid, int what);
        ET** getRemotAdressLocalSortOwnerOnly(size_t rid, int what);
        ET** getRemotAdressRemotSortOwnerOnly(size_t rid, int what);
        */
};

} // end namespace
#endif
