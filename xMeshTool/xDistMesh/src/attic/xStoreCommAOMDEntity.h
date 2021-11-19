/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef XSTOREAOMDENTITY_H
#define XSTOREAOMDENTITY_H

#include <set>

#include "xStoreCommEntity.h"

#include "mAttachableDataContainer.h"

namespace AOMD
{
    class mEntity;
}


namespace distmesh
{

/// Specialisation of xStoreCommEntity template class for AOMD::mEntity type 
//! It use AOMD attached data to store xStoreCommEntity
template <>
class xStoreCommEntity <AOMD::mEntity>
{
    public:
        /// Default constructor to initialise storage properly
        xStoreCommEntity(void);
        /// Copy constructor to initialize storage properly from an other instance
        xStoreCommEntity(const xStoreCommEntity& other);
        /// Destructor to clean corectelly storage
        ~xStoreCommEntity(void);
        /// Method to add in the storage information relate to a duplicate entity in a remote process
        //! \param [in] local adress of an entity duplicated in a other process
        //! \param [in] rproc remote pid 
        //! \param [in] radress remote adress of the duplicate entity
        void set(AOMD::mEntity *ladress, size_t rproc, AOMD::mEntity *radress);
        /// Method to get xCommEntity information related to entity e
        //!  \return xCommEntity relate to e
        //!  \return nul if no xCommEntity is associated to e
        xCommEntity * get(AOMD::mEntity* e);
        /// Remove xCommEntity information related to e in the container
        void clear(AOMD::mEntity* e);
        /// Remove  xCommConector information for rproc from xCommEntity related to e in the container
        //! If rproc was the last xCommConector in xCommEntity container then xCommEntity is removed
        //! from storage
        void clear(AOMD::mEntity* e, size_t rproc);
        /// Assign a proc_id information to instance
        //! \param [in] proc_id_ id of the current process this instance belong to.
        void assignProcId(size_t proc_id_);
        /// Tells if e entity is the owner among all proc of the information associated to it
        //! This is used to create communication pattern where only owner are sending to all
        //! related proc. This is removing intermediate communication in between non owner proc
        //! \param [in] e entity address which is supposed to have information to communicate
        //! \return true if e is the owner (This include entity with no associated xCommEntity)
        //! \return false if e is not the owner and have an associated xCommEntity
        bool isOwner(AOMD::mEntity*e);
    private:
        // Private class
        // Attachable data of type xCommEntity
        class xAttachableCommEntity : public AOMD::mAttachableData
        {
            public:
                ~xAttachableCommEntity();
                xCommEntity ace;
        };
        // Private methode
        // function to get attached connection
        xCommEntity * getAttachedCommEntity(AOMD::mEntity* e) const;
        /// function to attache a new connection to an entity removing old information
        void attachCommEntity(AOMD::mEntity* e, xCommEntity  & v);
        /// function to delete a connection from an entity 
        bool deleteCommEntity(AOMD::mEntity* e);
        // Private menber data
        unsigned int tag_bnd;
        size_t proc_id;
        bool proc_id_is_assigned;
        std::set<AOMD::mEntity*> attached;

};

} // end namspace
#endif
