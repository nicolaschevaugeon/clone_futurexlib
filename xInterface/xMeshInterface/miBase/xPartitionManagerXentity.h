#ifndef _XPARTITIONMANAGERXENTITY_H
#define _XPARTITIONMANAGERXENTITY_H

#include "xEntity.h"
#include "xEntityHashKey.h"
#include "xGeneralUnorderedMapDataManager.h"
#include "xPartitionManager.h"

namespace xinterface
{
namespace xmeshinterface
{
// cas particulier pour xEntity : la map est gérée par xEntity, mais on stocke uniquement une void* correspondant à la
// UniqueAdDress de xEntity. Ceci se fait avec la fonction insertUntyped(...) du partitionObject (po).
template <class DATATYPE>
using DataManagerXentity = xtool::xGeneralUnorderedMapDataManager<xEntity, DATATYPE, xEntityHashKey, xEntityEqualKey>;

//--------------------------------------------------------------------------
// classe de base
//--------------------------------------------------------------------------
class xPartitionManagerXentityBase : public xtool::xPartitionManager<DataManagerXentity>
{
  public:
   xPartitionManagerXentityBase(MPI_Comm _comm = MPI_COMM_WORLD) : xPartitionManager(_comm){};
   ~xPartitionManagerXentityBase() = default;
   virtual xEntity getXentityFromUniqueAddressAndQueryTag(const void* address, unsigned int tag) const = 0;
   virtual void print() const = 0;
   inline MPI_Comm getComm() const { return comm; };
};

class xPartitionManagerXentityUnion;

//--------------------------------------------------------------------------
// classe d'un partition_manager construit à partir d'une query
//--------------------------------------------------------------------------
class xPartitionManagerXentity : public xPartitionManagerXentityBase
{
  public:
   xPartitionManagerXentity(const xMeshQueryInterface& _query, MPI_Comm _comm = MPI_COMM_WORLD)
       : xPartitionManagerXentityBase(_comm), query(_query){};
   ~xPartitionManagerXentity() = default;
   xEntity getXentityFromUniqueAddressAndQueryTag(const void* address, unsigned int tag) const;
   inline const xMeshQueryInterface& getQuery() const { return query; };
   void print() const;

  protected:
   // case of simple partition_manager (not union), query is needed but
   xEntity getXentityFromUniqueAddress(const void* address) const;
   friend xPartitionManagerXentityUnion;

  private:
   const xMeshQueryInterface& query;
};

//----------------------------------------------------------------------------------
// classe d'un partition_manager construit à partir d'une union de plusieurs query :
//----------------------------------------------------------------------------------
// => utile pour un DOF manager s'appuyant sur  plusieurs queries/maillages
// xPartitionManagerXentityUnion----------------------------
class xPartitionManagerXentityUnion : public xPartitionManagerXentityBase
{
  private:
   std::vector<const xPartitionManagerXentity*> partman_container;

  public:
   xPartitionManagerXentityUnion(MPI_Comm _comm = MPI_COMM_WORLD) : xPartitionManagerXentityBase(_comm){};
   ~xPartitionManagerXentityUnion();
   void add(const xPartitionManagerXentity& new_element);
   const xtool::xConstPartitionObject<xEntity> getConstPartitionObject(xEntity& enti) const;
   xEntity getXentityFromUniqueAddressAndQueryTag(const void* address, unsigned int tag) const;
   void print() const;
};

}  // namespace xmeshinterface
}  // namespace xinterface

#endif
