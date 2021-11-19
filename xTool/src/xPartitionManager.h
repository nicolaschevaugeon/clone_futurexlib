/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _xPartitionManager_
#define _xPartitionManager_
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#include "mpi.h"
#include "xIteratorTools.h"

namespace xtool
{
class xUntypedRemoteObject;

template <class OBJECTTYPE>
class xRemoteObject;

/// xUntypedRemoteObjectsCollection represent a collection of xUntypedRemoteObject that describes copies of the same object on
/// different processors.
/** This class is internal to the library. User is not expected to  use it directly. \n
    An object can only have one copy on each of the processor of the group it belong too. This is why the operator < is needed for
   xUntypedRemoteObject \n The interface is stl like, we have the familliar, begin, end, empty functions.\n Rational and discution
   of the present implementation : \n More than one implementation of this concept is possible. We choosed here to implement it in
   terms of a sorted vector of  xUntypedRemoteObject. An alternative would be to store them on a set for example. The rational for
   the present choice is that we think that the memory foot print of the Collection should be as small as possible. Indeed, for
   the targetted applications, each mesh entity on the boundary beetwwen to process will store an instance of the
   xUntypedRemoteObjectsCollection. We might pay the price in terms of efficiency latter on : indeed the sorted vector approach
   migh involve too much of relocation and data mouvement for mesh entities on the boundary of a large numbre of partition. */
class xUntypedRemoteObjectsCollection
{
  public:
   typedef std::vector<xUntypedRemoteObject> container;
   typedef container::const_iterator const_iterator;
   typedef container::iterator iterator;
   typedef xRange<const_iterator> const_range;

   inline const_range range() const;
   inline const_iterator begin() const;
   inline const_iterator end() const;
   /// Return true if the remote collection is empty.
   inline bool empty() const;
   /// Return an iterator to the xUntypedRemoteObject on process rrank
   /** \param rrank : process id of the remote process for which we look for a xUntypedRemoteObject \n
       if no remote object on process rrank exist, find return an iterator to the end of the collection.  */
   inline const_iterator find(int rrank) const;
   /// Insert in the collection the address of the object on proc rrank.
   /**  if an object address already exist at for process rrank, it is replaced by the one given as argument.*/
   inline void insert(int rrank, const void *object_address);
   /// Remove from the collection the address of the object on proc rrank.
   /** if the address is not present in the Collection, the function does nothing.
       \param rrank : process id of the remote process from which one want to remove the remoteObject **/
   inline void remove(int rrank);
   inline iterator begin();
   inline iterator end();

  private:
   container remote_entities_container;
};

/// xConstPartitionObject is a class to represent an object that belong to the current partition and give access to the address of
/// it's remote copies if any.
/** Once a partition has been set up, this is the only object the user should really interact with. He is expected to get them
 from a xPartitionManager \n An object of this class contains a reference to the local object and ways to access to it's remote
 copies. \n a xConstPartitionObject has a unique owner. The "owner" is the process that "own" the object
 **/
template <class OBJECTTYPE>
class xConstPartitionObject
{
  public:
   typedef convertValueIterator<xUntypedRemoteObjectsCollection::const_iterator, xRemoteObject<OBJECTTYPE>> iter_t;
   typedef xRange<iter_t> const_range;
   /// isOwner return true if the object is owned by the current processor.
   /** If an object as no remote copies, it is owned by the current processor and the function return true */
   inline bool isOwner() const;
   inline bool hasRemoteObject() const { return !remote_collection.empty(); }
   /// getOwner return the xRemoteObject that own the object
   /** if the object as no remote copy or if the current processor own the object, getOwner return a xRemoteObject referring
       to the current proc and the local object. */
   inline xRemoteObject<OBJECTTYPE> getOwner() const;
   /// return a pointer to the address of the remote copy on process rrank.
   /**  if there is no remote copy on rrank, the function returns a nullptr \n
        if rrank == the id of the current process, getRemoteObjectOn return the address of the local object. */
   inline const OBJECTTYPE *getRemoteObjectOn(int rrank) const;
   inline const void *getRemoteUntypedObjectOn(int rrank) const;
   /// return the range of  remoteObject containing all the remotes copies of the object.
   /**  if there is no remote copy on rrank, the function returns an empty range */
   inline const_range getRemoteObjectsCollectionRange() const;
   template <template <class> class DATAMANAGER>
   friend class xPartitionManager;
   inline xConstPartitionObject(const OBJECTTYPE &_local_object, const xUntypedRemoteObjectsCollection &_remote_collection,
                                int _mpi_rank);

  private:
   /// xConstPartitionObject constructor is private : it is meant to be build by the xPartition manager
   // inline                      xConstPartitionObject (const OBJECTTYPE &_local_object,  const xUntypedRemoteObjectsCollection &
   // _remote_collection , int _mpi_rank);
   const OBJECTTYPE &local_object;
   const xUntypedRemoteObjectsCollection &remote_collection;
   int mpi_rank;
};

/// xPartitionObject is a class to represent an object that belong to the current partition and give access to the address of it's
/// remote copies if any.
/** it is similar to xConstPartitionObject, but also give the possibility to update the remote copie collection it gives acces to.
    An object of this class should only be use when setting up a partition. A user that just need to interogate a partition object
 in order to set up some communication should use a xConstPartitionObject instead. \n The user is expected to get it's
 xPartitionObject from a xPartitionManager. \n An object of this class contains a reference to the local object and ways to access
 to it's remote copies and remove or add remotes copies. \n a xPartitionObject has a unique owner. The "owner" is the process that
 "own" the object. \n an object of the xPartitionObject class always refer internally to a xPartitionManager. Any modification of
 the remoteObject collection modify the data contained by the xPartitionManager.
 **/
template <class OBJECTTYPE, class PARTITIONMANAGER>
class xPartitionObject
{
  public:
   typedef convertValueIterator<xUntypedRemoteObjectsCollection::const_iterator, xRemoteObject<OBJECTTYPE>> iter_t;
   typedef xRange<iter_t> const_range;
   /// See xConstPartitionObject
   inline bool isOwner() const;
   inline bool hasRemoteObject() const { return !(premote_collection->empty()); }
   /// See xConstPartitionObject
   inline xRemoteObject<OBJECTTYPE> getOwner() const;
   /// See xConstPartitionObject
   inline const OBJECTTYPE *getRemoteObjectOn(int rrank) const;
   /// See xConstPartitionObject
   inline const_range getRemoteObjectsCollectionRange() const;
   /// insert a new remote copy inside the remote object collection of the local object.
   /** \param rrank : process id of the process having the remote object.
       \param object_address a pointer to the address of the remote object. \n
       rrank must be in the [0, size[ i ntervalle, where size in the number of process in the in the MPI_Comm the xPartitionObject
      refer to. an assert check this in deubg mode. \n if rrank == current process, nothing is done. \n if a remote copy already
      exist in at pid rrank, it is replaced by this function. if the remoteObjectCollection was empty, a new one, associated to
      the local object is created in the xPartitionManager and the xPartitionObject now refer to it. **/
   inline void insert(int rrank, const OBJECTTYPE *object_address);
   /// remove the remote copy on process rrank from the remote object collection of the local object.
   inline void insertUntyped(int rrank, const void *object_address);
   /** it obviously don't delete the object copy on the remote processor (this would be non local) .
       This function is intended to be used when a remote processor as signaled that he just destroyed it's copy of the local
   object. \param rrank : process id of the process having the remote object. if rrank == current process, nothing is done. \n if
   no remote copy on process rrank is known by the remoteObjectCollection, nothing is done. \n if the remoteObjectCollection
   becomes empty, the data associated to the local object is removed from the xPartitionManager.
   **/
   inline void remove(int rrank);
   template <template <class> class DATAMANAGER>
   friend class xPartitionManager;

  private:
   /// private constructor. used only by the xPartitionManager
   inline xPartitionObject(OBJECTTYPE &_local_object, xUntypedRemoteObjectsCollection &_remote_collection,
                           PARTITIONMANAGER &_partman, int mpi_rank, int mpi_size);
   OBJECTTYPE &local_object;
   xUntypedRemoteObjectsCollection *premote_collection;
   PARTITIONMANAGER &partman;
   int mpi_rank;
   int mpi_size;
};

/// xPartitionManager instance are meant to store all needed datas to deal with object collection that might have remote copy
/// on different processor.
/** User interact with it mostly by calling the getConstPartitionObject and  getPartitionObject with a local object as an input
   parameter Internally it store xUntypedRemoteObjectsCollection with the help of a DATAMANAGER. A DATAMANAGER is a template class
    conforming to a simple interface to store data's associated to an object. \n
    exemples of DATAMANAGER conforming to the interface are given by the classes xinterface::aomd::xAttachedDataManagerAOMD and
   xUnorderedMapDataManager.\n The class xPartitionManager is template on a DATAMANAGER of xUntypedRemoteObjectsCollection, and
   privately derive fromit. From the user point of view, it is an implementation detail, but in fact it permits to customize a
   xPartitionManager to differents key<->data association mechanisms. */
template <template <class> class DATAMANAGER>
class xPartitionManager : private DATAMANAGER<xUntypedRemoteObjectsCollection>
{
  public:
   typedef xPartitionManager<DATAMANAGER> partman_t;
   typedef xUntypedRemoteObjectsCollection::const_iterator const_iterator;
   typedef xRange<const_iterator> const_range;
   typedef typename DATAMANAGER<xUntypedRemoteObjectsCollection>::c_iterKey_t c_iter_object_t;
   /// xPartitionManager constructor.
   /** \param _comm is the MPI_Comm refering to the processor group that will share the object refered to by the xPartition
      Manager. \n all the getConstPartitionObject and getPartitionObject will return xPartitionObject or xConstPartitionObject
      refering to the communicator. **/
   inline xPartitionManager(MPI_Comm _comm);
   /// From a const reference of OBJECTTYPE e, return a xConstPartitionObject that has e has it's local object.
   /** use this function to refer to the remotes copies of e **/
   template <class OBJECTTYPE>
   inline xConstPartitionObject<OBJECTTYPE> getConstPartitionObject(const OBJECTTYPE &e) const;
   /// From a reference of OBJECTTYPE e, return a xPartitionObject that has e has it's local object.
   /** use this function to refer to modify  the collections of  remotes copies of  e **/
   template <class OBJECTTYPE>
   inline xPartitionObject<OBJECTTYPE, partman_t> getPartitionObject(OBJECTTYPE &e);
   /// Cleanup all the data refering to e in the xPartitionManager. Use it to remove an object from the partition boundary
   template <class OBJECTTYPE>
   inline void remove(OBJECTTYPE &e);
   /// Fully cleanup the xPartitionManager. The only thing left is the MPI_Comm.
   /** This function actually call clear on the DATAMANAGER. **/
   inline void clear();
   /// From a const reference of OBJECTTYPE e, return true if any remote copie exist on other process for e.
   /** It return false if e have no remotes copies. It must be use only when you are not using remote information  whatever is the
      response. Otherwise getConstPartitionObject and getPartitionObject must be prefered has they give always an object. This
      object will tell you if there is remote copie or not. And if remote copies exist retriving information on them will be
      imediate (no query in the underliyng database of xPartitionManager) **/
   template <class OBJECTTYPE>
   inline bool hasRemoteCopy(const OBJECTTYPE &e) const;

   /// Give comunicator used with this partition manager
   inline MPI_Comm getComm() const;

   /// iterator on object having remote counterpart.
   inline c_iter_object_t beginObject() const;
   inline c_iter_object_t endObject() const;

   template <class OBJECTTYPE, class PARTITIONMANAGER>
   friend class xPartitionObject;

  private:
   template <class OBJECTTYPE>
   const xUntypedRemoteObjectsCollection *getRemoteObjectsCollection(const OBJECTTYPE &e) const;
   template <class OBJECTTYPE>
   xUntypedRemoteObjectsCollection *getRemoteObjectsCollection(OBJECTTYPE &e);
   template <class OBJECTTYPE>
   xUntypedRemoteObjectsCollection &setRemoteObjectsCollection(OBJECTTYPE &e);
   xUntypedRemoteObjectsCollection empty;

  protected:
   MPI_Comm comm;
   int mpi_rank;
   int mpi_size;
};

/// xUntypedRemoteObject is a class that describe an object of unknown type stored on a different process.
/**  This class is not meant to be used by library client. It's there to help the xPartitionManager do it's job.
    A xRemote object instance is meant to help refer to an objet that can exist on more than one process.\n
    The intent is to be able, with the help of an xPartitionMananger instance, to maintain and refer to copies
    of an object that live in more than one processor.
    It contain Two data : \n
    proc_id : the id of the mpi process that own a copy of your object. \n
    object_adress : a pointer to the address of the copy of your object on process proc_id. \n
    One never want to derefence the address stored in object address !! -> this address is only there to help
    communicate between copies of an entity. */
class xUntypedRemoteObject
{
  public:
   inline xUntypedRemoteObject(int _proc_id, const void *object_address) : proc_id(_proc_id), object_address(object_address) {}
   inline int getProcId() const { return proc_id; }
   inline const void *getObjectAddress() const { return object_address; }

  private:
   int proc_id;
   const void *object_address;
};

/// xRemoteObject is a template on OBJECTTYPE to represent an object on a remote processor
/** an object of this class store the process number and the address of an OBJECTTYPE on a remote processor.
    the address itself is never meant to be derefernced !  it's an address on a remote processor.
    It is used only to send message to it.
**/
template <class OBJECTTYPE>
class xRemoteObject
{
  public:
   /// xRemoteObject constructor.
   /** \param _proc_id represent the id of a remote processor,
       \param _object_address is a pointer storing the remote address of the object on proc proc_id
   **/
   inline xRemoteObject(int _proc_id, const OBJECTTYPE *_object_address);
   /// xRemoteObject constructor, from a xUntypedRemoteObject
   /** \param in : the input xUntypedRemoteObject \n
        xUntypedRemoteObject represent the same concept as xRemoteObject but the type of the remote object is unknown. \n
       (the address is stored as a void *). User of the library will probably never use this constructor. It is used inside
       the xPartitionObject and xConstPartitionObject iterators so that the user never has to manually cast from void to it's own
   type
   **/
   inline xRemoteObject(const xUntypedRemoteObject &in);
   /// return the process id on which the object live.
   inline int getProcId() const;
   /// return the address of the object in the remote process
   inline const OBJECTTYPE *getObjectAddress() const;
   inline const void *getUntypedObjectAddress() const;

  private:
   int proc_id;
   const OBJECTTYPE *object_address;
};

/// Comparaison operator between 2 xUntypedRemoteObject
/** This is not meant to be used by the client. This is really internal to the library.
    The intent is to be able to sort xUntypedRemoteObject describing copies of the same object on different processor.*/
inline bool operator<(const xUntypedRemoteObject &a, const xUntypedRemoteObject &b) { return a.getProcId() < b.getProcId(); }

/// amongs the local object pointed to by the range object_begin to object_end, count owmany belong to the current processor.
template <class ITER, template <class xUntypedRemoteObjectsCollection> class DATAMANAGER>
inline int countOwned(ITER object_begin, ITER object_end, const xPartitionManager<DATAMANAGER> &part_man);

/// check if the partitions boundary are correct. This is collective among all the processor of the MPI_Comm of the
/// xPartitionManager part_man
/** for each object in the list, look if it as remote entities. if it the case, check that these entities exist with the right
   address in remote processor
    return true if correct. return false if not. Write to std::cout message describing what is wrong if it's the case. **/
template <class ITERPENTITY, template <class xUntypedRemoteObjectsCollection> class DATAMANAGER>
bool checkPartition(ITERPENTITY itb, ITERPENTITY ite, const xPartitionManager<DATAMANAGER> &part_man);

}  // namespace xtool

namespace xtool
{
//***********************************************************
//***********************************************************
//***********************************************************
//     All Implementation start here !
//***********************************************************

//***********************************************************
//     xRemoteObject Implementation
//***********************************************************
template <class OBJECTTYPE>
inline xRemoteObject<OBJECTTYPE>::xRemoteObject(int _proc_id, const OBJECTTYPE *_object_address)
    : proc_id(_proc_id), object_address(_object_address)
{
}
template <class OBJECTTYPE>
inline xRemoteObject<OBJECTTYPE>::xRemoteObject(const xUntypedRemoteObject &in)
    : proc_id(in.getProcId()), object_address(static_cast<const OBJECTTYPE *>(in.getObjectAddress()))
{
}
template <class OBJECTTYPE>
inline int xRemoteObject<OBJECTTYPE>::getProcId() const
{
   return proc_id;
}
template <class OBJECTTYPE>
inline const OBJECTTYPE *xRemoteObject<OBJECTTYPE>::getObjectAddress() const
{
   return object_address;
}
// for xEntity, object_address is the key given as void*, not xEntity*
template <class OBJECTTYPE>
inline const void *xRemoteObject<OBJECTTYPE>::getUntypedObjectAddress() const
{
   return static_cast<const void *>(object_address);
}

//***********************************************************
//  xUntypedRemoteObjectsCollection Implementation
//***********************************************************
inline auto xUntypedRemoteObjectsCollection::range() const -> const_range { return const_range(begin(), end()); }

inline auto xUntypedRemoteObjectsCollection::begin() const -> const_iterator { return remote_entities_container.begin(); }
inline auto xUntypedRemoteObjectsCollection::end() const -> const_iterator { return remote_entities_container.end(); }
inline auto xUntypedRemoteObjectsCollection::begin() -> iterator { return remote_entities_container.begin(); }
inline auto xUntypedRemoteObjectsCollection::end() -> iterator { return remote_entities_container.end(); }
inline bool xUntypedRemoteObjectsCollection::empty() const { return remote_entities_container.empty(); }
inline auto xUntypedRemoteObjectsCollection::find(int rrank) const -> const_iterator
{
   auto it =
       std::lower_bound(remote_entities_container.begin(), remote_entities_container.end(), xUntypedRemoteObject(rrank, nullptr));
   if (it != remote_entities_container.end())
      if (it->getProcId() == rrank) return it;
   return remote_entities_container.end();
}

inline void xUntypedRemoteObjectsCollection::insert(int rrank, const void *object_address)
{
   xUntypedRemoteObject re(rrank, object_address);
   auto it = std::lower_bound(remote_entities_container.begin(), remote_entities_container.end(), re);
   if (it != remote_entities_container.end())
   {
      if (it->getProcId() == rrank) it = remote_entities_container.erase(it);
   }
   remote_entities_container.insert(it, re);
}

inline void xUntypedRemoteObjectsCollection::remove(int rrank)
{
   xUntypedRemoteObject re(rrank, nullptr);
   auto it = std::lower_bound(remote_entities_container.begin(), remote_entities_container.end(), re);
   if (it != remote_entities_container.end())
   {
      if (it->getProcId() == rrank) remote_entities_container.erase(it);
   }
}

//***********************************************************
// xConstPartitionObject Implementation
//***********************************************************
template <class OBJECTTYPE>
xConstPartitionObject<OBJECTTYPE>::xConstPartitionObject(const OBJECTTYPE &_local_object,
                                                         const xUntypedRemoteObjectsCollection &_remote_collection, int _mpi_rank)
    : local_object(_local_object), remote_collection(_remote_collection), mpi_rank(_mpi_rank)
{
}

template <class OBJECTTYPE>
inline bool xConstPartitionObject<OBJECTTYPE>::isOwner() const
{
   if (remote_collection.empty()) return true;
   return (remote_collection.begin()->getProcId() > mpi_rank);
};

template <class OBJECTTYPE>
inline xRemoteObject<OBJECTTYPE> xConstPartitionObject<OBJECTTYPE>::getOwner() const
{
   if (!remote_collection.empty())
   {
      xRemoteObject<OBJECTTYPE> ret(*(remote_collection.begin()));
      if (ret.getProcId() < mpi_rank) return ret;
   }
   return xRemoteObject<OBJECTTYPE>(mpi_rank, &local_object);
}

template <class OBJECTTYPE>
inline const OBJECTTYPE *xConstPartitionObject<OBJECTTYPE>::getRemoteObjectOn(int rrank) const
{
   if (mpi_rank == rrank) return &local_object;
   auto it = remote_collection.find(rrank);
   if (it != remote_collection.end()) return static_cast<const OBJECTTYPE *>((it->getObjectAddress()));
   return nullptr;
}
// especially for xEntity
template <class OBJECTTYPE>
inline const void *xConstPartitionObject<OBJECTTYPE>::getRemoteUntypedObjectOn(int rrank) const
{
   if (mpi_rank == rrank) return &local_object;
   auto it = remote_collection.find(rrank);
   if (it != remote_collection.end()) return static_cast<const void *>((it->getObjectAddress()));
   return nullptr;
}

template <class OBJECTTYPE>
inline auto xConstPartitionObject<OBJECTTYPE>::getRemoteObjectsCollectionRange() const -> const_range
{
   return const_range(iter_t(remote_collection.begin()), iter_t(remote_collection.end()));
};

//***********************************************************
// xPartitionObject Implementation
//***********************************************************
template <class OBJECTTYPE, class PARTITIONMANAGER>
xPartitionObject<OBJECTTYPE, PARTITIONMANAGER>::xPartitionObject(OBJECTTYPE &_local_object,
                                                                 xUntypedRemoteObjectsCollection &_remote_collection,
                                                                 PARTITIONMANAGER &_partman, int _mpi_rank, int _mpi_size)
    : local_object(_local_object),
      premote_collection(&_remote_collection),
      partman(_partman),
      mpi_rank(_mpi_rank),
      mpi_size(_mpi_size)
{
}
template <class OBJECTTYPE, class PARTITIONMANAGER>
inline bool xPartitionObject<OBJECTTYPE, PARTITIONMANAGER>::isOwner() const
{
   if (premote_collection->empty()) return true;
   return (premote_collection->begin()->getProcId() > mpi_rank);
}

template <class OBJECTTYPE, class PARTITIONMANAGER>
inline xRemoteObject<OBJECTTYPE> xPartitionObject<OBJECTTYPE, PARTITIONMANAGER>::getOwner() const
{
   if ((!premote_collection->empty()))
   {
      xRemoteObject<OBJECTTYPE> ret(*(premote_collection->begin()));
      if (ret.getProcId() < mpi_rank) return ret;
   }
   return xRemoteObject<OBJECTTYPE>(mpi_rank, &local_object);
}

template <class OBJECTTYPE, class PARTITIONMANAGER>
inline const OBJECTTYPE *xPartitionObject<OBJECTTYPE, PARTITIONMANAGER>::getRemoteObjectOn(int rrank) const
{
   if (mpi_rank == rrank) return &local_object;
   if (premote_collection->empty()) return nullptr;
   auto it = premote_collection->find(rrank);
   if (it != premote_collection->end()) return static_cast<const OBJECTTYPE *>((it->getObjectAddress()));
   return nullptr;
}

template <class OBJECTTYPE, class PARTITIONMANAGER>
inline auto xPartitionObject<OBJECTTYPE, PARTITIONMANAGER>::getRemoteObjectsCollectionRange() const -> const_range
{
   return const_range(iter_t(premote_collection->begin()), iter_t(premote_collection->end()));
};

template <class OBJECTTYPE, class PARTITIONMANAGER>
inline void xPartitionObject<OBJECTTYPE, PARTITIONMANAGER>::insert(int rrank, const OBJECTTYPE *object_address)
{
   if (rrank == mpi_rank) return;
   assert((rrank >= 0) && (rrank < mpi_size));
   if (premote_collection->empty())
   {
      premote_collection = &(partman.setRemoteObjectsCollection(local_object));
   }
   premote_collection->insert(rrank, object_address);
}

template <class OBJECTTYPE, class PARTITIONMANAGER>
inline void xPartitionObject<OBJECTTYPE, PARTITIONMANAGER>::remove(int rrank)
{
   if (premote_collection->empty()) return;
   if (rrank == mpi_rank) return;
   premote_collection->remove(rrank);
   if (premote_collection->empty())
   {
      partman.remove(local_object);
      premote_collection = &(partman.empty);
   }
}

//***********************************************************
// xPartitionManager Implementation
//***********************************************************
template <template <class> class DATAMANAGER>
inline xPartitionManager<DATAMANAGER>::xPartitionManager(MPI_Comm _comm) : comm(_comm)
{
   MPI_Comm_size(comm, &mpi_size);
   MPI_Comm_rank(comm, &mpi_rank);
}

template <template <class> class DATAMANAGER>
template <class OBJECTTYPE>
inline xConstPartitionObject<OBJECTTYPE> xPartitionManager<DATAMANAGER>::getConstPartitionObject(const OBJECTTYPE &e) const
{
   const xUntypedRemoteObjectsCollection *p_remobjcoll = getRemoteObjectsCollection(e);
   if (!p_remobjcoll)
   {
      return xConstPartitionObject<OBJECTTYPE>(e, empty, mpi_rank);
   }
   return xConstPartitionObject<OBJECTTYPE>(e, *p_remobjcoll, mpi_rank);
}

template <template <class> class DATAMANAGER>
template <class OBJECTTYPE>
xPartitionObject<OBJECTTYPE, xPartitionManager<DATAMANAGER>> xPartitionManager<DATAMANAGER>::getPartitionObject(OBJECTTYPE &e)
{
   xUntypedRemoteObjectsCollection *p_remobjcoll = getRemoteObjectsCollection(e);
   if (!p_remobjcoll)
   {
      return xPartitionObject<OBJECTTYPE, xPartitionManager<DATAMANAGER>>(e, empty, *this, mpi_rank, mpi_size);
   }
   return xPartitionObject<OBJECTTYPE, xPartitionManager<DATAMANAGER>>(e, *p_remobjcoll, *this, mpi_rank, mpi_size);
}

template <template <class> class DATAMANAGER>
template <class OBJECTTYPE>
inline void xPartitionManager<DATAMANAGER>::remove(OBJECTTYPE &e)
{
   DATAMANAGER<xUntypedRemoteObjectsCollection>::deleteData(e);
}

template <template <class> class DATAMANAGER>
inline void xPartitionManager<DATAMANAGER>::clear()
{
   DATAMANAGER<xUntypedRemoteObjectsCollection>::clear();
}

template <template <class> class DATAMANAGER>
template <class OBJECTTYPE>
inline bool xPartitionManager<DATAMANAGER>::hasRemoteCopy(const OBJECTTYPE &e) const
{
   return (getRemoteObjectsCollection(e) != nullptr);
}

template <template <class> class DATAMANAGER>
inline MPI_Comm xPartitionManager<DATAMANAGER>::getComm() const
{
   return comm;
}

template <template <class> class DATAMANAGER>
template <class OBJECTTYPE>
inline const xUntypedRemoteObjectsCollection *xPartitionManager<DATAMANAGER>::getRemoteObjectsCollection(
    const OBJECTTYPE &e) const
{
   return DATAMANAGER<xUntypedRemoteObjectsCollection>::getData(e);
}

template <template <class> class DATAMANAGER>
template <class OBJECTTYPE>
inline xUntypedRemoteObjectsCollection *xPartitionManager<DATAMANAGER>::getRemoteObjectsCollection(OBJECTTYPE &e)
{
   return DATAMANAGER<xUntypedRemoteObjectsCollection>::getData(e);
}

template <template <class> class DATAMANAGER>
template <class OBJECTTYPE>
inline xUntypedRemoteObjectsCollection &xPartitionManager<DATAMANAGER>::setRemoteObjectsCollection(OBJECTTYPE &e)
{
   return DATAMANAGER<xUntypedRemoteObjectsCollection>::setData(e);
}

template <template <class> class DATAMANAGER>
inline auto xPartitionManager<DATAMANAGER>::beginObject() const -> c_iter_object_t
{
   return this->beginKey();
}
template <template <class> class DATAMANAGER>
inline auto xPartitionManager<DATAMANAGER>::endObject() const -> c_iter_object_t
{
   return this->endKey();
}
//***********************************************************
// countOwned  Implementation
//***********************************************************
template <class ITER, template <class xUntypedRemoteObjectsCollection> class DATAMANAGER>
inline int countOwned(ITER entities_begin, ITER entities_end, const xPartitionManager<DATAMANAGER> &part_man)
{
   int nowned = 0;
   for (auto it = entities_begin; it != entities_end; ++it)
   {
      auto partobj = part_man.getConstPartitionObject(*(*it));
      if (partobj.isOwner()) ++nowned;
   }
   return nowned;
}

//***********************************************************
// checkPartition Implementation. uses the following sparse_sendreceive for the implementation of the passing message.
//            sparse_sendreceive should be replaced by more generics function of the library when available.
//***********************************************************
template <class EXCHANGEDATA>
void sparse_sendreceive(const std::map<int, std::vector<EXCHANGEDATA>> &sendbuf,
                        std::map<int, std::vector<EXCHANGEDATA>> &recvbuf, int commtag, MPI_Comm comm)
{
   // const bool debug = false;
   int mpisize, rank;
   MPI_Comm_size(comm, &mpisize);
   MPI_Comm_rank(comm, &rank);

   // send size of data to receive
   std::vector<int> receive_sizes_sendbuf(mpisize, 0);
   for (auto i : sendbuf) receive_sizes_sendbuf[i.first] = i.second.size();
   std::vector<int> receive_sizes_receivebuf(mpisize, 0);
   MPI_Alltoall(receive_sizes_sendbuf.data(), 1, MPI_INT, receive_sizes_receivebuf.data(), 1, MPI_INT, comm);

   for (int i = 0; i < mpisize; ++i)
   {
      if (rank == i)
      {
         for (auto sendto : sendbuf)
         {
            int to = sendto.first;
            if (to != i)
               MPI_Send(sendto.second.data(), sendto.second.size() * sizeof(EXCHANGEDATA), MPI_BYTE, to, commtag, comm);
            else
               recvbuf[to] = sendto.second;
         }
      }
      else
      {
         if (receive_sizes_receivebuf[i] != 0)
         {
            recvbuf[i].resize(receive_sizes_receivebuf[i]);
            MPI_Status stat;
            MPI_Recv(recvbuf[i].data(), receive_sizes_receivebuf[i] * sizeof(EXCHANGEDATA), MPI_BYTE, i, commtag, comm, &stat);
         }
      }
   }
   MPI_Barrier(comm);
}

template <class ITERPENTITY, template <class xUntypedRemoteObjectsCollection> class DATAMANAGER>
bool checkPartition(ITERPENTITY itb, ITERPENTITY ite, const xPartitionManager<DATAMANAGER> &part_man)
{
   int rank, mpisize;
   MPI_Comm_size(part_man.getComm(), &mpisize);
   MPI_Comm_rank(part_man.getComm(), &rank);

   typedef typename ITERPENTITY::value_type pentity;
   typedef typename std::remove_pointer<pentity>::type entity;
   struct data
   {
      const entity *lv;
      const entity *rv;
   };
   std::map<int, std::vector<data>> send;
   std::map<int, std::vector<data>> recv;
   for (auto it = itb; it != ite; ++it)
   {
      entity *e = *it;
      for (auto r : part_man.getConstPartitionObject(*e).getRemoteObjectsCollectionRange())
      {
         if (r.getObjectAddress() == nullptr)
         {
            std::cout << "P" << rank << " Entity " << e << " has null remote " << std::endl;
            return false;
         }
         data datae = {e, r.getObjectAddress()};
         send[r.getProcId()].push_back(datae);
      }
   }

   sparse_sendreceive(send, recv, 11225, part_man.getComm());

   for (auto recv_from : recv)
   {
      int from = recv_from.first;
      for (auto data : recv_from.second)
      {
         const entity *le = data.rv;
         const entity *re = data.lv;
         const entity *re_inpm = part_man.getConstPartitionObject(*le).getRemoteObjectOn(from);
         if (re != re_inpm)
         {
            std::cout << "P" << rank << " Receive from P" << from << " that entity "
                      << "P" << from << " " << re << " is "
                      << "P" << rank << " " << le << " While P" << rank << " Think it's P" << from << " " << re_inpm << std::endl;
            return false;
         }
      }
   }
   return true;
}

}  // namespace xtool

#endif
