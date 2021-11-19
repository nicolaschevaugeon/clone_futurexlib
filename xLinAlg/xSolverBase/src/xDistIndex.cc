/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
// std
#include "xDistIndex.h"

#include <iostream>
#include <limits>
#include <utility>

namespace xlinalg
{
//===== initIndexKeysInfoManagers ============================================================================================
class initIndexKeysInfoManagers
{
  public:
   // mandatory types for information class familly and key manager class familly
   typedef xtool::homogeneous_data_style_trait data_style_trait;
   typedef xtool::send_only_keys_communication_trait communication_trait;
   typedef const xDistIndex::idx_t *information_key_t;
   typedef std::pair<information_key_t, xDistIndex::idx_t> information_t;

   /// Constructor
   initIndexKeysInfoManagers(xDistIndex &dist_index_) : dist_index(dist_index_) {}

   // mandatory method for information class familly
   information_key_t localObjectKey(const xDistIndex::idx_t &o) { return dist_index.getPackedAdress(o); }
   information_t getInfo(information_key_t key, int sendto)
   {
      // the key is the information
      return std::make_pair(key, *key);
   }
   void setInfo(const std::vector<information_t> &info, int receivedfrom)
   {
      for (auto inf : info)
      {
         xDistIndex::idx_t *p = const_cast<xDistIndex::idx_t *>(dist_index.getPackedAdress(inf.second));
         auto po = dist_index.id_part_man.getPartitionObject(*p);
         po.insert(receivedfrom, inf.first);
      }
   }
   std::set<int> getMessageRanks(const information_key_t &lk) { return dist_index.to_from[lk - dist_index.packed_id.data()]; }
   information_t one_instance;

  private:
   xDistIndex &dist_index;
};

//===== xDistIndex ============================================================================================
xDistIndex::xDistIndex(MPI_Comm world_)
    : id_part_man(world_),
      world(world_),
      nb_local_idx(0),
      nb_glob_idx(0),
      min_idx(std::numeric_limits<idx_t>::max()),
      max_idx(0),
      p_global_id(nullptr),
      to_be_finalized(true),
      inc_ordering(false)
{
}
xDistIndex::xDistIndex(const xDistIndex &other)
    : packed_id(other.packed_id),
      global_id(other.global_id),
      all_nb_local_idx(other.all_nb_local_idx),
      all_disp_local_idx(other.all_disp_local_idx),
      id_part_man(other.id_part_man),
      world(other.world),
      nb_local_idx(other.nb_local_idx),
      nb_glob_idx(other.nb_glob_idx),
      min_idx(other.min_idx),
      max_idx(other.max_idx),
      p_global_id(&global_id[0] - min_idx),
      to_be_finalized(false),
      inc_ordering(other.inc_ordering)
{
   if (other.to_be_finalized)
   {
      std::cout << "Ouch ! you can't copy a xDistIndex instance it is not finalize ! In " << __FILE__ << " " << __LINE__ << " "
                << __DATE__ << " " << __TIME__ << std::endl;
      throw -1;
   }
}
void xDistIndex::insertIndex(idx_t idx)
{
   assert(idx > 0);
   // check finalized status
   assert(to_be_finalized);

   // in global
   if (idx < min_idx) min_idx = idx;
   if (static_cast<idx_t>(global_id.size()) < idx)
   {
      global_id.resize(idx, -1);
      p_global_id = &global_id[0] - 1;
   }
   p_global_id[idx] = packed_id.size();

   // in packed
   packed_id.push_back(idx);
}
void xDistIndex::insertToFrom(idx_t idx, int tf)
{
   assert(idx > 0);
   // check finalized status
   assert(to_be_finalized);

   to_from.resize(packed_id.size());
   assert(p_global_id[idx] < packed_id.size());
   assert(p_global_id[idx] > -1);
   to_from[p_global_id[idx]].insert(tf);
}
void xDistIndex::finalize(idx_t nb_local, idx_t max_local_idx, bool inc_ordering_)
{
   int proc_id, nb_proc;
   // getting rank in world
   MPI_Comm_rank(world, &proc_id);
   // getting size of world
   MPI_Comm_size(world, &nb_proc);

   // setting nb local
   assert(!(nb_local > packed_id.size()));
   nb_local_idx = nb_local;
   all_nb_local_idx.resize(nb_proc);
   MPI_Allgather(&nb_local_idx, sizeof(idx_t), MPI_BYTE, all_nb_local_idx.data(), sizeof(idx_t), MPI_BYTE, world);

   // setting disp local
   all_disp_local_idx.resize(nb_proc + 1);
   int k = 1;
   all_disp_local_idx[0] = 0;
   for (auto n : all_nb_local_idx)
   {
      all_disp_local_idx[k] = all_disp_local_idx[k - 1] + n;
      ++k;
   }

   // set numdof global to all proc
   max_idx = max_local_idx;
   MPI_Allreduce((void *)&max_local_idx, (void *)&nb_glob_idx, 1, MPI_INT, MPI_MAX, world);

   // finalize indexes
   idx_t mxg = global_id.size();
   if (mxg)
   {
      if (min_idx > 1)
      {
         std::copy(global_id.begin() + min_idx - 1, global_id.end(), global_id.begin());
         mxg = mxg - min_idx + 1;
         global_id.resize(mxg);
      }
      p_global_id = &global_id[0] - min_idx;
   }
   else
      p_global_id = global_id.data();

   if (nb_proc > 1)
   {
      idx_t mx = packed_id.size();
      if (mx != static_cast<idx_t>(to_from.size()))
      {
#ifndef NDEBUG
         std::cout << "In " << __FILE__ << " " << __LINE__ << " " << __DATE__ << " " << __TIME__ << std::endl;
         std::cout
             << "Ouch ! you try to finalize index and partition manager but size of to_from is not eguale to size of packed_id."
             << std::endl;
         std::cout << "Hop it is not the mess ! resize" << std::endl;
#endif
         to_from.resize(packed_id.size());
      }
      initIndexKeysInfoManagers key_info_manager(*this);
      xtool::xKeyContainerSendOrRecv<initIndexKeysInfoManagers::information_key_t> key_container(world);
      // We need for exchanger to say it is finalized (use of getPackedIndex)
      // We use unprotected function during this loop but we are sure that those call are safe because they
      // are donne on entered values.
      to_be_finalized = false;
      for (idx_t i = 0; i < mxg; ++i)
      {
         if (global_id[i] > -1)
         {
            key_container.accumulateKeys(packed_id[global_id[i]], key_info_manager);  // accumulate keys to be exchanged
         }
      }
      exchangeInformation(key_container, key_info_manager);
      to_from.clear();
   }
   else
      to_be_finalized = false;

   // checking of increasing order
   inc_ordering = inc_ordering_;
   if (inc_ordering)
   {
      // TODO : with nb_local,packed_id,all_disp_local_idx and min_idx things can be done like checking
      // that there is no holes and index start in conformance with disp. This may help changing inc_ordering to false.
      // But as inc_ordering==false is not done for now we don't anything here also.
      ;  //
   }

   // if not in increasing order create helping information
   if (!inc_ordering)
   {
      ;  // TODO create somthing that help skiping holes and folow index in increasing order
   }
}
xDistIndex::idx_t xDistIndex::getGlobalIndexSize() const
{
   assert(!to_be_finalized);
   return nb_glob_idx;
}
xDistIndex::idx_t xDistIndex::getPackedIndexSize() const
{
   assert(!to_be_finalized);
   return packed_id.size();
}
xDistIndex::idx_t xDistIndex::getPackedLocalIndexSize() const
{
   assert(!to_be_finalized);
   return nb_local_idx;
}
const xDistIndex::idx_t *xDistIndex::getAllPackedLocalIndexSize() const
{
   assert(!to_be_finalized);
   return all_nb_local_idx.data();
}
const xDistIndex::idx_t *xDistIndex::getAllPackedLocalIndexDisp() const
{
   assert(!to_be_finalized);
   return all_disp_local_idx.data();
}
xDistIndex::idx_t xDistIndex::getPackedIndex(const idx_t idx) const
{
   assert(!to_be_finalized);
   assert(!(idx < min_idx));
   idx_t i = p_global_id[idx];
   assert(i > -1);
   return i;
}
bool xDistIndex::testIndexPresent(idx_t idx) const
{
   assert(!to_be_finalized);
   if (idx < min_idx || idx > max_idx || p_global_id[idx] < 0) return false;
   return true;
}
const xDistIndex::idx_t *xDistIndex::getPackedAdress(idx_t idx) const { return &packed_id[getPackedIndex(idx)]; }

auto xDistIndex::begin() const -> idx_container_const_iterator_t { return packed_id.begin(); }
auto xDistIndex::end() const -> idx_container_const_iterator_t { return packed_id.end(); }
MPI_Comm xDistIndex::getComm() const { return world; }
xtool::xConstPartitionObject<xDistIndex::idx_t> xDistIndex::getConstPartitionObject(const xDistIndex::idx_t &idx) const
{
   const idx_t *p = getPackedAdress(idx);
   return id_part_man.getConstPartitionObject(*p);
}
bool xDistIndex::isInIncreasingIndexing() const { return inc_ordering; }
}  // namespace xlinalg
