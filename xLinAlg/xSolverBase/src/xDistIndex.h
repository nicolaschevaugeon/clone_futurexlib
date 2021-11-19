/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _XDISTINDEX_H
#define _XDISTINDEX_H

// std
#include <set>
#include <vector>

// parallel
#include "mpi.h"

// eXlibris_tools
#include "xDataExchanger.h"
#include "xPartitionManager.h"
#include "xUnorderedMapDataManager.h"

namespace xlinalg
{
class initIndexKeysInfoManagers;

/*!
 * This class store index in a packed distributed format. It offers method to store indexes and retrieve there packed position.
 */
class xDistIndex
{
  public:
   /*! This is the type of index used in this class.
    *
    * One would say, hey ! why not using size_t !
    * In fact we are using this class in the context of MPI, and in this norm its "int" which is used.
    * And when people want more they use compiler to transform int4 in int8.
    * So a compatibility choice.
    *
    * This choice offer nice guard : negative value. This is simplifying implementation from original version
    */
   typedef int idx_t;

  private:
   typedef std::vector<idx_t> idx_container_t;
   template <class T>
   using data_manager_t = xtool::xUnorderedMapDataManager<idx_container_t::value_type, T>;

  public:
   typedef xtool::xPartitionManager<data_manager_t> partman_t;
   typedef std::vector<idx_t>::const_iterator idx_container_const_iterator_t;

   /// Unique explicit constructor
   //! \param [in] world_ communicator on which distributed index have been constructed
   xDistIndex(MPI_Comm world_);

   /*! Default copy constructor is modified to :
    *  - avoid copy of a non finalized index
    *  - set correctly internal address
    *  \param [in] other an instance of an other xDistIndex.
    */
   xDistIndex(const xDistIndex &other);

   /*! This method add an index in the storage.
    * The order you call this method is important. Says if you insert index 5 before index 1 this
    * will says that in packed form index 5 is before index 1.
    * The first entered index are supposed to be owned by this proc (see finalize ). Here owned refer to owned concept
    * of xPartitionManager class. Mainly an index will be owned by the proc having the smallest id in the communicator
    * \param [in] idx index added to storage (Fortran indexing)
    */
   void insertIndex(idx_t idx);

   /*! This method set a connection for a  given index in the storage.
    * This method have to be called on a index only if this index have been stored. A checking is done in debug version.
    * Index may be added after last call to this method.
    * \param [in] idx index already in the storage   (Fortran indexing)
    * \param [in] tf process id in world_ which have also this index
    */
   void insertToFrom(idx_t idx, int tf);

   /*! This method finalize object creation. After calling this method you will not be able to store any
    *  other index (assert are put in insert method). Some communication (heavy?) are done in this fucntion to set internal
    * feature.
    * \param [in] nb_local is the number of owned index by this proc. This should correspond to the first nb_local call
    * to insertIndex.
    * \param [in] max_local_idx is the greatest index value entered so far for this proc.
    * \param [in] inc_ordering_ tells if indexing is in increasing order of proc without holes or not. For now there is nothing
    * behind this parameter. It is used thought by xDistVector class. It is set to false by default, as false mean real general
    * indexing. There can be holes and indexing may jump back in some proc ... The idea is that this setting is done by the class
    * in the future.
    */
   void finalize(idx_t nb_local, idx_t max_local_idx, bool inc_ordering_ = false);

   /*! This method give the total number of index that may have been stored in all proc.
    *  This size may be greater then the real number of index really stored because it correspond to max_local_idx. User
    *  may have some index not stored, i.e. holes in there indexation.
    */
   idx_t getGlobalIndexSize() const;

   /*! This method give the total number of index entered in this instance (proc).
    */
   idx_t getPackedIndexSize() const;

   /*! This method give the total number of index entered in this instance which are local/owned by this proc.
    * This correspond to nb_local argument of finalize method
    */
   idx_t getPackedLocalIndexSize() const;

   /*! This method give the total number of index entered in this instance which are local/owned by all  proc.
    * This correspond to nb_local argument of finalize method given in all proc
    * \return A pointer to an array of local size. Array has the size of the communicator.
    */
   const idx_t *getAllPackedLocalIndexSize() const;

   /*! Having a fictuous global vector grouping all local index (in packed form) of all proc (in increasing proc order),
    * this method give all the displacement from begin of this vector for each proc.
    * \return A pointer to an array of displacement. Array has the size of the communicator.
    */
   const idx_t *getAllPackedLocalIndexDisp() const;

   /*! This method give the communicator on which index are set
    */
   MPI_Comm getComm() const;

   /*! This method give the beginning of the packed index container
    * \return A  const iterator of type idx_container_const_iterator_t
    */
   idx_container_const_iterator_t begin() const;
   /*! This method give the past end of the packed index container
    * \return A const  iterator of type idx_container_const_iterator_t
    */
   idx_container_const_iterator_t end() const;

   /*! This method give the index of the idx parameter in the packed index container
    * \param [in] idx index (in fortran indexing) for which we look for it's location in packed index container
    * \return A index/location in  packed index container
    */
   idx_t getPackedIndex(const idx_t idx) const;

   /*! This method give the address of the idx parameter in the packed index container
    * \param [in] idx index (in fortran indexing) for which we look for it's location in packed index container
    * \return A address in  packed index container
    */
   const idx_t *getPackedAdress(idx_t idx) const;

   /*! This method test if the index idx parameter is stored in instance at local level
    * (i.e. in the packed index container )
    *  It return true if found and false otherwise.
    * \param [in] idx index (in fortran indexing) for which we test for presence
    * \return A bool which is true if idx is found
    */
   bool testIndexPresent(idx_t idx) const;

   /*! This method give the const partition object related to index  idx
    * \param [in] idx index (in Fortran indexing) for which we look for it's associated xConstPartitionObject
    * \return A xConstPartitionObject related to index parameter
    */
   xtool::xConstPartitionObject<xDistIndex::idx_t> getConstPartitionObject(const xDistIndex::idx_t &idx) const;

   /*! This method return true if indexing is in increasing order on increasing proc id with no holes
    * \return A boolean which is true if in increasing order
    * Note: For now this is set by user in finalize method.
    */
   bool isInIncreasingIndexing() const;

  private:
   friend class initIndexKeysInfoManagers;
   idx_container_t packed_id;
   idx_container_t global_id;
   idx_container_t all_nb_local_idx;
   idx_container_t all_disp_local_idx;
   std::vector<std::set<int>> to_from;
   partman_t id_part_man;
   MPI_Comm world;
   idx_t nb_local_idx, nb_glob_idx;
   idx_t min_idx;
   idx_t max_idx;
   idx_t *p_global_id;
   bool to_be_finalized;
   bool inc_ordering;
};

}  // namespace xlinalg

#endif
