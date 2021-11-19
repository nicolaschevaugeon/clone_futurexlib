/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
// standard includes
#include <algorithm>
#include <fstream>
#include <functional>
#include <vector>
// Trellis/AOMD
#include "mEdge.h"
#include "mFace.h"
#include "mMesh.h"
#include "mRegion.h"
#include "mTet.h"
#include "mVertex.h"
#include "modeler.h"
// Trellis/model
#include "GEntity.h"
// Accumulate
#include <numeric>
// Macro
#include "xMacro.h"

namespace xmeshtool
{
inline auto range(const AOMD::mMesh &m, int dim) -> xtool::xRange<decltype(m.begin(dim))>
{
   assert(dim >= 0 && dim <= 3);
   return xtool::xRange<decltype(m.begin(dim))>(m.begin(dim), m.end(dim));
}

class exchangerKeyManagerSimpleSendAndReceiveKeysFromPartitionManagerAOMD
{
  public:
   exchangerKeyManagerSimpleSendAndReceiveKeysFromPartitionManagerAOMD(const partitionManager &_partman) : partman(_partman) {}
   typedef const AOMD::mEntity *information_key_t;
   information_key_t localObjectKey(const AOMD::mEntity &o) const { return &o; }
   information_key_t remoteObjectKey(const xtool::xRemoteObject<AOMD::mEntity> &ro, const AOMD::mEntity &lo) const
   {
      return ro.getObjectAddress();
   }
   xtool::xConstPartitionObject<AOMD::mEntity> getConstPartitionObject(const AOMD::mEntity &e) const
   {
      return partman.getConstPartitionObject(e);
   }

  protected:
   const partitionManager &partman;
};

// strategy :  all proc first count the vertices they own, then exchanges this numbers with an MPI_Scan. Now all the process know
// where to start the numbering of there owned vertices. numbering is done. then not owned vertices retrieve there number from
// owner proc during exchange.
//
template <class ITER>
inline xinterface::aomd::xAttachedDataManagerAOMD<int> globalEntitiesNumbering(const ITER &b, const ITER &e,
                                                                               const partitionManager &part_man)
{
   xtool::xRange<ITER> entities_range(b, e);
   const int nbowned = countOwned(b, e, part_man);
   const MPI_Comm comm = part_man.getComm();

   int iend;
   MPI_Scan(const_cast<int *>(&nbowned), &iend, 1, MPI_INT, MPI_SUM, comm);
   const int istart = iend - nbowned;

   xinterface::aomd::xAttachedDataManagerAOMD<int> entities_numbering;

   int i = istart;
   for (auto p_entity : entities_range)
   {
      if (part_man.getConstPartitionObject(*p_entity).isOwner()) entities_numbering.setData(*p_entity) = i++;
   }

   class exchangerInfoManagerNumberEntities
   {
     public:
      exchangerInfoManagerNumberEntities(xinterface::aomd::xAttachedDataManagerAOMD<int> &_entities_numbering)
          : entities_numbering(_entities_numbering)
      {
      }
      typedef xtool::homogeneous_data_style_trait data_style_trait EXLIBRIS_MACRO_WARNUNUSEDTYPE;
      typedef xtool::send_and_recv_keys_communication_trait communication_trait EXLIBRIS_MACRO_WARNUNUSEDTYPE;
      typedef const AOMD::mEntity *information_key_t;
      typedef int information_t;

      information_t getInfo(information_key_t key, int sendto)
      {
         const AOMD::mEntity &e = *key;
         const int *i = entities_numbering.getData(e);
         assert(i);
         return *i;
      };

      void setInfo(information_key_t key, const information_t &i, int receivedfrom)
      {
         AOMD::mEntity &e = const_cast<AOMD::mEntity &>(*key);
         assert(!entities_numbering.getData(e));  // entity e should not have yet a number at this point.
         entities_numbering.setData(e) = i;
      };

     private:
      xinterface::aomd::xAttachedDataManagerAOMD<int> &entities_numbering;
   };

   typedef const AOMD::mEntity *key_t;
   typedef exchangerInfoManagerNumberEntities exchangerInfoManager;
   typedef exchangerKeyManagerSimpleSendAndReceiveKeysFromPartitionManagerAOMD exchangerKeyManager;
   typedef xtool::xKeyContainerSendAndRecv<key_t> exchangerKeyContainer;

   exchangerInfoManager exchanger_info_manager(entities_numbering);
   exchangerKeyManager exchanger_key_manager(part_man);
   exchangerKeyContainer exchanger_key_container(part_man.getComm());

   exchanger_key_container.accumulateKeysOwnerScatter(b, e, exchanger_key_manager);
   exchangeInformation(exchanger_key_container, exchanger_info_manager);
   return entities_numbering;
}

inline xinterface::aomd::xAttachedDataManagerAOMD<int> globalVerticesNumbering(AOMD::mMesh &m, const partitionManager &part_man)
{
   auto vrange = range(m, 0);
   return globalEntitiesNumbering(vrange.begin(), vrange.end(), part_man);
}

inline bool isTopLevel(const AOMD::mEntity &e)
{
   const int level = e.getLevel();
   if (level == 3)
      return true;
   else
      return !(e.size(level + 1));
}

template <int NBW>
xinterface::aomd::xAttachedDataManagerAOMD<int> setParmetisPartitionWithWeight(
    AOMD::mMesh &m, const partitionManager &part_man, int nb_part,
    xinterface::aomd::xAttachedDataManagerAOMD<std::array<xinterface::parmetis::ParMetisInterface::parmetis_indx_t, NBW>>
        &weights)
{
   xinterface::aomd::xAttachedDataManagerAOMD<int> global_vertex_id = globalVerticesNumbering(m, part_man);
   std::array<int, 4> nb_top_dim = {{0, 0, 0, 0}};
   for (int dim = 0; dim < 3; ++dim)
   {
      for (auto pentity : range(m, dim))
      {
         if (isTopLevel(*pentity)) ++nb_top_dim[dim];
      }
   }
   nb_top_dim[3] = m.size(3);
   const int nb_top = std::accumulate(nb_top_dim.begin(), nb_top_dim.end(), 0);
   MPI_Comm comm = part_man.getComm();
   int rank;
   MPI_Comm_rank(comm, &rank);
   int color = (nb_top) ? 1 : 0;
   MPI_Comm newcomm;
   MPI_Comm_split(comm, color, rank, &newcomm);
   xinterface::aomd::xAttachedDataManagerAOMD<int> futurelem_rank;
   if (color)
   {
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t ncon = NBW;
      int new_mpi_size;
      MPI_Comm_size(newcomm, &new_mpi_size);
      int iend;

      MPI_Scan(const_cast<int *>(&nb_top), &iend, 1, MPI_INT, MPI_SUM, newcomm);
      int istart = iend;
      std::vector<xinterface::parmetis::ParMetisInterface::parmetis_indx_t> elmdist(new_mpi_size + 1, 0);
      MPI_Allgather(&istart, 1, MPI_INT, elmdist.data() + 1, 1, MPI_INT, newcomm);

      std::vector<xinterface::parmetis::ParMetisInterface::parmetis_indx_t> eptr(nb_top + 1);
      std::vector<xinterface::parmetis::ParMetisInterface::parmetis_indx_t> eind(nb_top_dim[0] + 2 * nb_top_dim[1] +
                                                                                 3 * nb_top_dim[2] + 4 * nb_top_dim[3]);
      std::vector<xinterface::parmetis::ParMetisInterface::parmetis_indx_t> elmwgt(nb_top * ncon);
      eptr[0] = 0;
      int pos = 1;
      int weight_idx = 0;
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t index = 0;
      for (int dim = 0; dim <= 3; ++dim)
      {
         for (auto pentity : range(m, dim))
         {
            if (isTopLevel(*pentity))
            {
               if (dim == 0)
               {
                  const int *pnodeid = global_vertex_id.getData(*pentity);
                  assert(pnodeid);
                  eind[index++] = *pnodeid;
               }
               else
               {
                  for (int in = 0; in < dim + 1; ++in)
                  {
                     const int *pnodeid = global_vertex_id.getData(*(pentity->get(0, in)));
                     assert(pnodeid);
                     eind[index++] = *pnodeid;
                  }
               }
               eptr[pos++] = index;
               auto *pw = weights.getData(*pentity);
               xinterface::parmetis::ParMetisInterface::parmetis_indx_t *pew = &elmwgt[weight_idx];
               if (pw)
               {
                  auto &rw = *pw;
                  for (int i = 0; i < NBW; ++i) pew[i] = rw[i];
               }
               else
               {
                  cout << "Entity ...TODO... got no weight => null by default" << endl;
                  for (int i = 0; i < NBW; ++i) pew[i] = 0;
               }
               weight_idx += NBW;
            }
         }
      }

      xinterface::parmetis::ParMetisInterface::parmetis_indx_t wgtflag = 2;
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t numflag = 0;
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t ncommonnodes = 3;
      if (nb_top_dim[2]) ncommonnodes = 2;
      if (nb_top_dim[1]) ncommonnodes = 1;
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t nparts = nb_part;

      std::vector<xinterface::parmetis::ParMetisInterface::parmetis_real_t> tpwgts(ncon * nparts, 1.f / nparts);
      std::vector<xinterface::parmetis::ParMetisInterface::parmetis_real_t> ubvec(ncon, 1.05);
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t edgecut;

      std::vector<xinterface::parmetis::ParMetisInterface::parmetis_indx_t> part(nb_top, -1);
      std::array<xinterface::parmetis::ParMetisInterface::parmetis_indx_t, 3> options = {{1, 0, 0}};
      //    std::cout<< "start Parmetis " << std::endl;
      xinterface::parmetis::ParMetisInterface::PartMeshKway(elmdist.data(), eptr.data(), eind.data(), elmwgt.data(), wgtflag,
                                                            numflag, ncon, ncommonnodes, nparts, tpwgts.data(), ubvec.data(),
                                                            options.data(), newcomm, edgecut, part.data());
      int i = 0;
      for (int dim = 0; dim <= 3; ++dim)
      {
         for (auto pentity : range(m, dim))
         {
            if (isTopLevel(*pentity))
            {
               futurelem_rank.setData(*pentity) = part[i++];
            }
         }
      }
   }
   //  std::cout<< "PARMETIS Done "<< "Edge cut "<< edgecut << std::endl;
   return futurelem_rank;
}

}  // namespace xmeshtool
