/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __VALUE__MANAGER__H
#define __VALUE__MANAGER__H

#include "xDataExchanger.h"
#include "xDataExchangerTools.h"
#include "xPartitionManager.h"
#include "xValKey.h"
#include "xValKeyDataManager.h"
#include "xValManager.h"
#include "xValue.h"

namespace xfem
{
template <typename VT>
class xValueManagerDist : public xValManager<xValKey, xValue<VT>, xHashValKey, xEqualValKey>
{
  private:
   using mesh_partman_t = partmanAOMD_t;
   template <typename T>
   using data_manager_t = xValKeyDataManager<T>;
   //        template < typename T >
   //        using data_manager_t = xtool::xGeneralUnorderedMapDataManager < xValKey, T,xHashValKey, xEqualValKey >;
  public:
   typedef xtool::xPartitionManager<data_manager_t> partman_t;
   typedef VT value_t;
   typedef xValue<VT> xvalue_t;

   /// Constructor that store a reference to the  manager for which we
   //! will need distributed information
   xValueManagerDist();

   /// Destructor needed  to clen dynamic allocation
   ~xValueManagerDist();

   /// function that create distributed information
   void genPartitionManager(const mesh_partman_t &mesh_partman);

   /// function to obtain the const version of the distributed information
   auto getPartitionManager(void) const -> const partman_t &;
   /// function to obtain the non const version of the distributed information
   auto getPartitionManager(void) -> partman_t &;

  private:
   partman_t *partition_manager;

   // =======================================================================================================================================
   // Local class to exchange information
   class xDistKeyKeymanager
   {
     public:
      // mandatory traits
      typedef xValKey information_key_t;

      xDistKeyKeymanager(const mesh_partman_t &part_man_) : part_man(part_man_) {}

      // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
      //        typedef std::iterator_traits < xValueManagerDist<VT>::map_iterator >::value_type data_t;
      using data_t = typename std::iterator_traits<typename xValueManagerDist<VT>::map_iterator>::value_type;

      // mandatory methods
      information_key_t localObjectKey(const data_t &o) { return o.first; }
      std::set<int> getMessageRanks(const information_key_t &key)
      {
         // key is the xValkey to treat

         // first get entity of this key => local object
         AOMD::mEntity *lo = key.getEnti();

         // get remote object
         xtool::xConstPartitionObject<AOMD::mEntity> po = part_man.getConstPartitionObject(*lo);

         // int container (will be returned by this method), empty by default (no comunication)
         std::set<int> res;

         // if object have remote copy
         if (po.hasRemoteObject())
         {
            // grabe remote proc id for object and store them in res
            for (const auto ro : po.getRemoteObjectsCollectionRange())
            {
               res.insert(ro.getProcId());
            }
         }

         return res;
      }

     private:
      const mesh_partman_t &part_man;
   };
   // this structure is here to transfert the xvalkey information with its origine
   struct transXvalKey
   {
      xValKey *adress;
      AOMD::mEntity *recever_Enti;
      int Refe;
      xValKey::ids_size_t Phys;
      xValKey::ids_size_t Geom;
   };

   class xDistKeyInfoManager
   {
     public:
      typedef xtool::homogeneous_data_style_trait data_style_trait;
      typedef xtool::send_only_keys_communication_trait communication_trait;

      using information_key_t = typename xDistKeyKeymanager::information_key_t;
      typedef transXvalKey information_t;

      xDistKeyInfoManager(xValueManagerDist<VT> &dm_, const mesh_partman_t &mesh_partman_)
          : dm(dm_), part_man(dm.getPartitionManager()), mesh_part_man(mesh_partman_)
      {
      }

      // mandatory methods
      information_t getInfo(const information_key_t &key, int sendto)
      {
         // Get adresse in dm of the key
         const xValKey *lka = dm.findKeyAdress(key);
         assert(lka);

         // first get entity of this key to create remote
         AOMD::mEntity *lo = key.getEnti();

         // get remote object
         xtool::xConstPartitionObject<AOMD::mEntity> po = mesh_part_man.getConstPartitionObject(*lo);

         assert(po.hasRemoteObject());
         transXvalKey tmp;
         tmp.adress = const_cast<xValKey *>(lka);
         tmp.recever_Enti = const_cast<AOMD::mEntity *>(po.getRemoteObjectOn(sendto));
         tmp.Refe = key.getRefe();
         tmp.Phys = key.getPhys();
         tmp.Geom = key.getGeom();
         return tmp;
      }
      void setInfo(const std::vector<information_t> &infos, int receivedfrom)
      {
         for (auto info : infos)
         {
            xValKey tmp(info.Phys, info.Geom, info.recever_Enti, info.Refe);

            // if this key is present in dm then a connection must be set
            if (dm.findKeyAdress(tmp))
            {
               auto npo = part_man.getPartitionObject(tmp);
               npo.insert(receivedfrom, info.adress);
            }
         }
      }

     private:
      xValueManagerDist<VT> &dm;
      typename xValueManagerDist<VT>::partman_t &part_man;
      const mesh_partman_t &mesh_part_man;
   };
};

// using xDoubleManager = xValueManagerDist<double>;
// using xFloatManager = xValueManagerDist<float>;
// using xDoubleComplexManager = xValueManagerDist<std::complex<double> >;

// =======================================================================================================================================
// xValueManagerDist implementation

template <typename VT>
xValueManagerDist<VT>::xValueManagerDist() : partition_manager(nullptr)
{
}

template <typename VT>
xValueManagerDist<VT>::~xValueManagerDist()
{
   if (partition_manager) delete partition_manager;
   partition_manager = nullptr;
}

/// function that create distributed informations
template <typename VT>
void xValueManagerDist<VT>::genPartitionManager(const mesh_partman_t &mesh_partman)
{
   MPI_Comm world = mesh_partman.getComm();
   if (partition_manager)
      partition_manager->clear();
   else
      partition_manager = new partman_t(world);

   // set keys: xValkey that are on some proc boundary
   xDistKeyKeymanager key_manager(mesh_partman);
   xtool::xKeyContainerSendOrRecv<typename xDistKeyKeymanager::information_key_t, xHashValKey, xEqualValKey> key_container(world);
   key_container.accumulateKeys(this->MapContainer.begin(), this->MapContainer.end(), key_manager);

   // exchange: creation of partition_manager
   xDistKeyInfoManager info_manager(*this, mesh_partman);
   xtool::exchangeInformation(key_container, info_manager);
}

template <typename VT>
auto xValueManagerDist<VT>::getPartitionManager(void) const -> const partman_t &
{
   assert(partition_manager);
   return *partition_manager;
}

template <typename VT>
auto xValueManagerDist<VT>::getPartitionManager(void) -> partman_t &
{
   assert(partition_manager);
   return *partition_manager;
}

template <typename VT>
class xValueLinearCombination;

template <typename VT>
class xCloneTemplateValue
{
  public:
   xValue<VT> *operator()(const xValue<VT> *val, const std::map<xValue<VT> *, xValue<VT> *> &coresp) const
   {
      xValue<VT> *new_val = nullptr;

      if (typeid(*val) == typeid(xSingleValue<VT>))
         new_val = xCloneValue(static_cast<const xSingleValue<VT> *>(val),
                               (static_cast<const std::map<xValue<VT> *, xValue<VT> *> &>(coresp)));
      else if (typeid(*val) == typeid(xValueLinearCombination<VT>))
         new_val = xCloneValue(static_cast<const xValueLinearCombination<VT> *>(val),
                               (static_cast<const std::map<xValue<VT> *, xValue<VT> *> &>(coresp)));
      else
      {
         std::cout << "Value of type " << typeid(*val).name() << " is not (yet) clonable !" << std::endl;
         throw -12345;
      }

      return new_val;
   }
};

}  // namespace xfem

#endif
