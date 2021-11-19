/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
// DistMesh
#include "xLoadBalanceTools.h"

namespace xmeshtool
{
class ranksData
{
  public:
   std::vector<int> getNewRanks() const
   {
      std::vector<int> new_ranks = futur_ranks;
      auto new_end =
          std::set_difference(futur_ranks.begin(), futur_ranks.end(), old_ranks.begin(), old_ranks.end(), new_ranks.begin());
      new_ranks.resize(new_end - new_ranks.begin());
      return new_ranks;
   }
   std::vector<int> getRemoveRanks() const
   {
      std::vector<int> rm_ranks = old_ranks;
      auto rm_end =
          std::set_difference(old_ranks.begin(), old_ranks.end(), futur_ranks.begin(), futur_ranks.end(), rm_ranks.begin());
      rm_ranks.resize(rm_end - rm_ranks.begin());
      return rm_ranks;
   }

   const std::vector<int> &getOldRanks() const { return old_ranks; }
   const std::vector<int> &getFuturRanks() const { return futur_ranks; }
   void addToOldRanks(int pid)
   {
      auto it = std::lower_bound(old_ranks.begin(), old_ranks.end(), pid);
      if (it == old_ranks.end())
      {
         old_ranks.insert(it, pid);
         return;
      }
      if (pid != *it) old_ranks.insert(it, pid);
   }
   void addToFuturRanks(int pid)
   {
      auto it = std::lower_bound(futur_ranks.begin(), futur_ranks.end(), pid);
      if (it == futur_ranks.end())
      {
         futur_ranks.insert(it, pid);
         return;
      }
      if (pid != *it) futur_ranks.insert(it, pid);
   }

  private:
   std::vector<int> old_ranks;
   std::vector<int> futur_ranks;
};

typedef xinterface::aomd::xAttachedDataManagerAOMD<ranksData> ranksDataManager;

class exchangerInfoManagerEntitiesRanks
{
  public:
   exchangerInfoManagerEntitiesRanks(ranksDataManager &_rank_manager) : rank_manager(_rank_manager){};
   typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
   typedef xtool::send_and_recv_keys_communication_trait communication_trait;
   typedef const AOMD::mEntity *information_key_t;

   void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto)
   {
      //  std::cout << "getRanks Info" << std::endl;
      const AOMD::mEntity &e = *key;
      // std::cout << " Entity " << &e << std::endl;
      assert(rank_manager.getData(e));
      const auto &rank_data = (*rank_manager.getData(e));
      const auto &futur_ranks = rank_data.getFuturRanks();
      const int nb = futur_ranks.size();
      buff.pack(&nb, 1, MPI_INT);
      buff.pack(futur_ranks.data(), nb, MPI_INT);
      rank_manager.deleteData(const_cast<AOMD::mEntity &>(e));
   }

   void setInfo(information_key_t key, const xtool::xMpiOutputBuffer &buff, int receivedfrom)
   {
      // std::cout << "setRanks Info" << std::endl;
      int nb;
      buff.unPack(&nb, 1, MPI_INT);
      std::vector<int> rankscontiguous(nb);
      buff.unPack(rankscontiguous.data(), nb, MPI_INT);
      AOMD::mEntity &e = const_cast<AOMD::mEntity &>(*key);
      auto &rank_data = rank_manager.setData(e);
      for (auto pid : rankscontiguous) rank_data.addToFuturRanks(pid);
   };

   size_t getApproxDataSize() const { return 3 * sizeof(int); }

  private:
   ranksDataManager &rank_manager;
};

ranksDataManager setRankManager(const AOMD::mMesh &m, const partitionManager &part_man,
                                const xinterface::aomd::xAttachedDataManagerAOMD<int> &targetproc)
{
   int rank;
   MPI_Comm_rank(part_man.getComm(), &rank);
   ranksDataManager ranks_manager;
   xtool::xKeyContainerSendAndRecv<const AOMD::mEntity *> exchanger_key_container(part_man.getComm());
   exchangerKeyManagerSimpleSendAndReceiveKeysFromPartitionManagerAOMD exchanger_key_manager(part_man);
   exchangerInfoManagerEntitiesRanks exchanger_info_manager_ranks(ranks_manager);
   const std::array<std::string, 4> ltostring = {{"vertex", "edge", "face", "region"}};
   for (int dim = 3; dim >= 0; --dim)
   {
      for (auto it_pentity = m.begin(dim); it_pentity != m.end(dim); ++it_pentity)
      {
         AOMD::mEntity &entity = (*(*(it_pentity)));
         // std::cout << " Entity " << &entity << std::endl;
         auto po_entity = part_man.getConstPartitionObject(entity);
         auto &rank_data_entity = ranks_manager.setData(entity);
         if (po_entity.isOwner())
         {
            //	std::cout << " isOWNER " << std::endl;
            rank_data_entity.addToOldRanks(rank);
            for (auto ro_entity : po_entity.getRemoteObjectsCollectionRange())
            {
               rank_data_entity.addToOldRanks(ro_entity.getProcId());
            }
            if (isTopLevel(entity))
            {
               // std::cout << " istopLevel " << std::endl;
               const int *pfutur_rank = targetproc.getData(entity);
               if (!pfutur_rank)
               {
                  std::cout << "ERROR : " << ltostring[dim] << " " << &entity << " Does not have any data  in targetproc  "
                            << __FILE__ << ":" << __LINE__ << std::endl;
                  throw;
               }
               rank_data_entity.addToFuturRanks(*pfutur_rank);
            }
         }
         if (dim > 0)
         {
            // std::cout << "ENTITY " << &entity << std::endl;
            for (auto it_pdown_entity = entity.begin(dim - 1); it_pdown_entity != entity.end(dim - 1); ++it_pdown_entity)
            {
               auto &down_entity = *(*it_pdown_entity);
               // std::cout << " Down ENTITY " << &down_entity;
               auto &rank_data_down_entity = ranks_manager.setData(down_entity);
               for (auto futur_rank : rank_data_entity.getFuturRanks())
               {
                  // std::cout << " " << futur_rank << std::endl;
                  rank_data_down_entity.addToFuturRanks(futur_rank);
               }
               // std::cout<< std::endl;
            }
         }
      }

      if (dim > 0)
      {
         // std::cout << "CHECKING rank dim-1 "<< dim-1  << " : "<< std::endl;
         // for(auto it = m.begin(dim-1); it != m.end(dim-1); ++it){
         //	const AOMD::mEntity & e = *(*it);
         //	std::cout << "checking frank before exchange of " << &e << " : ";
         //	auto * rank_data_e = ranks_manager.getData(e);
         //	assert(rank_data_e);
         //	for (int r : rank_data_e->getFuturRanks()){
         //	  std::cout << " " << r;
         //	}
         //	std::cout << std::endl;
         exchanger_key_container.accumulateKeysOwnerGather(m.begin(dim - 1), m.end(dim - 1), exchanger_key_manager);
         exchangeInformation(exchanger_key_container, exchanger_info_manager_ranks);
         exchanger_key_container.clearKeys();
      }
   }
   // std::cout <<"## Will exit " << std::endl;
   return ranks_manager;
}

class exchangerKeyManagerSimpleSendOrReceiveKeysFromPartitionManagerAOMD
{
  public:
   typedef const AOMD::mEntity *information_key_t;
   information_key_t localObjectKey(const AOMD::mEntity &o) const { return &o; }
   virtual std::set<int> getMessageRanks(information_key_t lk) const = 0;
};

class exchangerKeyManagerCreateEntities : public exchangerKeyManagerSimpleSendOrReceiveKeysFromPartitionManagerAOMD
{
  public:
   exchangerKeyManagerCreateEntities(const partitionManager &_partman, const ranksDataManager &_rank_manager)
       : partman(_partman), rank_manager(_rank_manager)
   {
   }

   std::set<int> getMessageRanks(information_key_t lk) const override
   {
      std::set<int> targets;
      const AOMD::mEntity &entity = *lk;
      auto po_entity = partman.getConstPartitionObject(entity);
      if (po_entity.isOwner())
      {
         const auto *prank_data_entity = rank_manager.getData(entity);
         assert(prank_data_entity);
         const auto &ranks_data_entity = *prank_data_entity;
         auto ranks_to_rcopy = ranks_data_entity.getNewRanks();
         std::copy(ranks_to_rcopy.begin(), ranks_to_rcopy.end(), std::inserter(targets, targets.begin()));
      }
      return targets;
   }

  protected:
   const partitionManager &partman;
   const ranksDataManager &rank_manager;
};

class exchangerInfoManagerEntitiesBase
{
  protected:
   void reallyPackMoreDataEntity(const AOMD::mEntity &entity, xtool::xMpiInputBuffer &buff, int sendto) const
   {
      pmore_data_entity->getInfo(entity, buff, sendto);
      return;
   }
   void reallyUnPackMoreDataEntity(AOMD::mEntity &entity, const xtool::xMpiOutputBuffer &buff, int receivedfrom)
   {
      pmore_data_entity->setInfo(entity, buff, receivedfrom);
      return;
   }
   void reallyDelMoreDataEntity(AOMD::mEntity &entity)
   {
      pmore_data_entity->delEntity(entity);
      return;
   }
   void notReallyPackMoreDataEntity(const AOMD::mEntity &entity, xtool::xMpiInputBuffer &buff, int sendto) const { return; }
   void notReallyUnPackMoreDataEntity(AOMD::mEntity &entity, const xtool::xMpiOutputBuffer &buff, int receivedfrom) { return; }
   void notReallyDelMoreDataEntity(AOMD::mEntity &entity) { return; }

  public:
   exchangerInfoManagerEntitiesBase(AOMD::mMesh &_m, partitionManager &_part_man, ranksDataManager &_rank_manager,
                                    moreDataEntity *_pmore_data_entity)
       : m(_m), part_man(_part_man), rank_manager(_rank_manager), pmore_data_entity(_pmore_data_entity)
   {
      setMpiRoType();
      MPI_Comm_rank(part_man.getComm(), &rank);
      if (pmore_data_entity)
      {
         packMoreDataEntity = std::bind(&exchangerInfoManagerEntitiesBase::reallyPackMoreDataEntity, this, std::placeholders::_1,
                                        std::placeholders::_2, std::placeholders::_3);
         unPackMoreDataEntity = std::bind(&exchangerInfoManagerEntitiesBase::reallyUnPackMoreDataEntity, this,
                                          std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
         delMoreDataEntity = std::bind(&exchangerInfoManagerEntitiesBase::reallyDelMoreDataEntity, this, std::placeholders::_1);
      }
      else
      {
         packMoreDataEntity = std::bind(&exchangerInfoManagerEntitiesBase::notReallyPackMoreDataEntity, this,
                                        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
         unPackMoreDataEntity = std::bind(&exchangerInfoManagerEntitiesBase::notReallyUnPackMoreDataEntity, this,
                                          std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
         delMoreDataEntity =
             std::bind(&exchangerInfoManagerEntitiesBase::notReallyDelMoreDataEntity, this, std::placeholders::_1);
      }
   }

   ~exchangerInfoManagerEntitiesBase() { MPI_Type_free(&mpi_ro_t); }
   struct ro_t
   {
      int id;
      AOMD::mEntity *ra;
   };

   MPI_Datatype mpi_ro_t;
   void setMpiRoType()
   {
      ro_t test;
      int count = 2;
      int array_of_blocklengths[] = {1, 1};
      MPI_Aint a0;
      MPI_Aint aid;
      MPI_Aint ara;
      MPI_Get_address(&test, &a0);
      MPI_Get_address(&test.id, &aid);
      MPI_Get_address(&test.ra, &ara);
      MPI_Aint array_of_displacements[] = {aid - a0, ara - a0};
      MPI_Datatype array_of_types[] = {MPI_INT, MPI_AINT};
      MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &mpi_ro_t);
      MPI_Type_commit(&mpi_ro_t);
   }
   void packGEntity(const AOMD::GEntity &g, xtool::xMpiInputBuffer &buff) const
   {
      std::array<int, 2> dim_tag = {{g.dim(), g.tag()}};
      buff.pack(dim_tag.data(), 2, MPI_INT);
   }

   AOMD::GEntity *unPackGEntity(const xtool::xMpiOutputBuffer &buff) const
   {
      std::array<int, 2> dim_tag;
      buff.unPack(dim_tag.data(), 2, MPI_INT);
      return m.getGEntity(dim_tag[1], dim_tag[0]);
   }

   void packRemotesEntities(const AOMD::mEntity &entity, xtool::xMpiInputBuffer &buff) const
   {
      const auto po = part_man.getConstPartitionObject(entity);
      const auto rrange = po.getRemoteObjectsCollectionRange();
      const int nbr = rrange.size() + 1;
      buff.pack(&nbr, 1, MPI_INT);
      std::vector<ro_t> ro_buff(nbr);
      ro_buff[0] = {rank, const_cast<AOMD::mEntity *>(&entity)};
      int i = 1;
      for (auto ro : rrange) ro_buff[i++] = {ro.getProcId(), const_cast<AOMD::mEntity *>(ro.getObjectAddress())};
      buff.pack(ro_buff.data(), nbr, mpi_ro_t);
   }

   void unPackRemotesEntities(const xtool::xMpiOutputBuffer &buff, AOMD::mEntity &entity)
   {
      auto po = part_man.getPartitionObject(entity);
      int nbr;
      buff.unPack(&nbr, 1, MPI_INT);
      std::vector<ro_t> ro_buff(nbr);
      buff.unPack(ro_buff.data(), nbr, mpi_ro_t);
      for (auto ro : ro_buff) po.insert(ro.id, ro.ra);
   }

  public:
   std::function<void(const AOMD::mEntity &, xtool::xMpiInputBuffer &, int)> packMoreDataEntity;
   std::function<void(AOMD::mEntity &, const xtool::xMpiOutputBuffer &, int)> unPackMoreDataEntity;
   std::function<void(AOMD::mEntity &)> delMoreDataEntity;
   AOMD::mMesh &m;
   partitionManager &part_man;
   ranksDataManager &rank_manager;
   moreDataEntity *pmore_data_entity;
   std::map<int, std::vector<const AOMD::mEntity *>> created_entities_from;
   int rank;
};

class exchangerInfoManagerCreateEntities
{
  public:
   typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
   typedef xtool::send_only_keys_communication_trait communication_trait;
   typedef const AOMD::mEntity *information_key_t;
   exchangerInfoManagerCreateEntities(exchangerInfoManagerEntitiesBase &_base) : base(_base){};
   virtual void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto) const = 0;
   virtual void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom) = 0;
   virtual size_t getApproxDataSize() const = 0;

  protected:
   exchangerInfoManagerEntitiesBase &base;
};

class exchangerInfoManagerCreateVerticesBase : public exchangerInfoManagerCreateEntities
{
  public:
   typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
   typedef xtool::send_only_keys_communication_trait communication_trait;
   typedef const AOMD::mEntity *information_key_t;

   void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto) const override
   {
      const AOMD::mVertex &vertex = static_cast<const AOMD::mVertex &>(*key);
      Trellis_Util::mPoint p = vertex.point();
      buff.pack(&p(0), 3, MPI_DOUBLE);
      getVertexInfo(vertex, buff);
      base.packGEntity(*vertex.getClassification(), buff);
      base.packRemotesEntities(vertex, buff);
      base.packMoreDataEntity(vertex, buff, sendto);
   }
   void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom) override
   {
      while (!buff.exhausted())
      {
         AOMD::mVertex *pv = setVertex(buff);
         base.created_entities_from[receivedfrom].push_back(pv);
         base.unPackRemotesEntities(buff, *pv);
         base.unPackMoreDataEntity(*pv, buff, receivedfrom);
      }
   }
   size_t getApproxDataSize() const override { return approx_size; }

  protected:
   exchangerInfoManagerCreateVerticesBase(exchangerInfoManagerEntitiesBase &_base)
       : exchangerInfoManagerCreateEntities(_base), approx_size(0){};

   std::function<AOMD::mVertex *(const xtool::xMpiOutputBuffer &)> setVertex;
   std::function<void(const AOMD::mVertex &, xtool::xMpiInputBuffer &)> getVertexInfo;
   size_t approx_size;
};

class exchangerInfoManagerCreateVertices : public exchangerInfoManagerCreateVerticesBase
{
  public:
   exchangerInfoManagerCreateVertices(exchangerInfoManagerEntitiesBase &_base, bool fixed_vertex_id)
       : exchangerInfoManagerCreateVerticesBase(_base)
   {
      if (fixed_vertex_id)
      {
         setVertex = std::bind(&exchangerInfoManagerCreateVertices::setVertexFix, this, std::placeholders::_1);
         getVertexInfo =
             std::bind(&exchangerInfoManagerCreateVertices::getVertexInfoFix, this, std::placeholders::_1, std::placeholders::_2);
         // 3 double for X,Y, Z, 2 int for classification, 1 int for id , 1 int nbrentity nbrentity  * remote entity_size
         // rentitysiet = 1 int 1,
         approx_size = 3 * sizeof(double) + 2 * sizeof(int) + sizeof(int) + sizeof(int) +
                       1 * (sizeof(int) + sizeof(int *));  //+ pmore_data_entity->getApproxDataSize();
      }
      else
      {
         setVertex = std::bind(&exchangerInfoManagerCreateVertices::setVertexFree, this, std::placeholders::_1);
         getVertexInfo = std::bind(&exchangerInfoManagerCreateVertices::getVertexInfoFree, this, std::placeholders::_1,
                                   std::placeholders::_2);
         // 3 double for X,Y, Z, 2 int for classification , 1 int nbrentity nbrentity  * remote entity_size rentitysiet = 1 int 1,
         approx_size = 3 * sizeof(double) + 2 * sizeof(int) + sizeof(int) +
                       1 * (sizeof(int) + sizeof(int *));  //+ pmore_data_entity->getApproxDataSize();
      }
   };

   void getVertexInfoFree(const AOMD::mVertex &vertex, xtool::xMpiInputBuffer &buff) const { ; }
   void getVertexInfoFix(const AOMD::mVertex &vertex, xtool::xMpiInputBuffer &buff) const
   {
      int id = vertex.getId();
      buff.pack(&id, 1, MPI_INT);
   }
   AOMD::mVertex *setVertexFree(const xtool::xMpiOutputBuffer &buff)
   {
      std::array<double, 3> coord;
      buff.unPack(coord.data(), 3, MPI_DOUBLE);
      AOMD::GEntity *pg = base.unPackGEntity(buff);
      return base.m.createVertex(coord[0], coord[1], coord[2], pg);
   }
   AOMD::mVertex *setVertexFix(const xtool::xMpiOutputBuffer &buff)
   {
      std::array<double, 3> coord;
      int id;
      buff.unPack(coord.data(), 3, MPI_DOUBLE);
      buff.unPack(&id, 1, MPI_INT);
      AOMD::GEntity *pg = base.unPackGEntity(buff);
      return base.m.createVertex(id, coord[0], coord[1], coord[2], pg);
   }
};

class exchangerInfoManagerCreateEdges : public exchangerInfoManagerCreateEntities
{
  public:
   exchangerInfoManagerCreateEdges(exchangerInfoManagerEntitiesBase &_base) : exchangerInfoManagerCreateEntities(_base){};
   void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto) const override
   {
      const AOMD::mEdge &edge = static_cast<const AOMD::mEdge &>(*key);
      std::array<AOMD::mVertex *, 2> rov_array;
      for (int i = 0; i < 2; ++i)
      {
         const AOMD::mVertex &v = static_cast<const AOMD::mVertex &>(*edge.get(0, i));
         auto pov = base.part_man.getConstPartitionObject(v);
         rov_array[i] = const_cast<AOMD::mVertex *>(pov.getRemoteObjectOn(sendto));
         assert(rov_array[i]);
      }
      buff.pack(rov_array.data(), 2, MPI_AINT);
      base.packGEntity(*edge.getClassification(), buff);
      base.packRemotesEntities(edge, buff);
      base.packMoreDataEntity(edge, buff, sendto);
   }

   void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom) override
   {
      while (!buff.exhausted())
      {
         std::array<AOMD::mVertex *, 2> pvarray;
         buff.unPack(pvarray.data(), 2, MPI_AINT);
         AOMD::GEntity *pg = base.unPackGEntity(buff);
         AOMD::mEdge *pe = base.m.createEdge_oneLevel(pvarray[0], pvarray[1], pg);
         base.created_entities_from[receivedfrom].push_back(pe);
         base.unPackRemotesEntities(buff, *pe);
         base.unPackMoreDataEntity(*pe, buff, receivedfrom);
      }
   }
   size_t getApproxDataSize() const override { return sizeof(void *) * 2 + sizeof(int) * 2 + (sizeof(int) + sizeof(void *)) * 2; }
};

class exchangerInfoManagerCreateTriangles : public exchangerInfoManagerCreateEntities
{
  public:
   exchangerInfoManagerCreateTriangles(exchangerInfoManagerEntitiesBase &_base) : exchangerInfoManagerCreateEntities(_base){};
   void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto) const override
   {
      const AOMD::mFace &face = static_cast<const AOMD::mFace &>(*key);
      assert(face.size(0) == 3);
      std::array<AOMD::mEdge *, 3> roe_array;
      for (int i = 0; i < 3; ++i)
      {
         const AOMD::mEdge &e = static_cast<const AOMD::mEdge &>(*face.get(1, i));
         auto poe = base.part_man.getConstPartitionObject(e);
         roe_array[i] = const_cast<AOMD::mEdge *>(poe.getRemoteObjectOn(sendto));
         assert(roe_array[i]);
      }
      buff.pack(roe_array.data(), 3, MPI_AINT);
      base.packGEntity(*face.getClassification(), buff);
      base.packRemotesEntities(face, buff);
      base.packMoreDataEntity(face, buff, sendto);
   }

   void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom) override
   {
      while (!buff.exhausted())
      {
         std::array<AOMD::mEdge *, 3> fearray;
         buff.unPack(fearray.data(), 3, MPI_AINT);
         AOMD::GEntity *pg = base.unPackGEntity(buff);
         AOMD::mFace *pf = base.m.createFaceWithEdges_oneLevel(fearray[0], fearray[1], fearray[2], pg);
         base.created_entities_from[receivedfrom].push_back(pf);
         base.unPackRemotesEntities(buff, *pf);
         base.unPackMoreDataEntity(*pf, buff, receivedfrom);
      }
   }

   size_t getApproxDataSize() const override { return sizeof(void *) * 3 + sizeof(int) * 2 + (sizeof(int) + sizeof(void *)) * 2; }
};

class exchangerInfoManagerCreateTetrahedra : public exchangerInfoManagerCreateEntities
{
  public:
   exchangerInfoManagerCreateTetrahedra(exchangerInfoManagerEntitiesBase &_base) : exchangerInfoManagerCreateEntities(_base){};
   void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto) const override
   {
      const AOMD::mRegion &region = static_cast<const AOMD::mRegion &>(*key);
      assert(region.size(0) == 4);
      std::array<AOMD::mFace *, 4> rof_array;
      for (int i = 0; i < 4; ++i)
      {
         const AOMD::mFace &f = static_cast<const AOMD::mFace &>(*region.get(2, i));
         auto pov = base.part_man.getConstPartitionObject(f);
         rof_array[i] = const_cast<AOMD::mFace *>(pov.getRemoteObjectOn(sendto));
         assert(rof_array[i]);
      }
      buff.pack(rof_array.data(), 4, MPI_AINT);
      base.packGEntity(*region.getClassification(), buff);
      base.packRemotesEntities(region, buff);
      base.packMoreDataEntity(region, buff, sendto);
   }

   void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom) override
   {
      while (!buff.exhausted())
      {
         std::array<AOMD::mFace *, 4> pfarray;
         buff.unPack(pfarray.data(), 4, MPI_AINT);
         AOMD::GEntity *pg = base.unPackGEntity(buff);
         AOMD::mRegion *pr = base.m.createTetWithFaces_oneLevel(pfarray[0], pfarray[1], pfarray[2], pfarray[3], pg);
         base.created_entities_from[receivedfrom].push_back(pr);
         base.unPackRemotesEntities(buff, *pr);
         base.unPackMoreDataEntity(*pr, buff, receivedfrom);
      }
   }

   size_t getApproxDataSize() const override { return sizeof(void *) * 4 + sizeof(int) * 2 + (sizeof(int) + sizeof(void *)) * 2; }
};

class exchangerInfoManagerSendBackEntitiesAddresses
{
  public:
   typedef xtool::homogeneous_data_style_trait data_style_trait;
   typedef xtool::recv_only_keys_communication_trait communication_trait;
   typedef const AOMD::mEntity *information_key_t;
   typedef const AOMD::mEntity *information_t;

   exchangerInfoManagerSendBackEntitiesAddresses(exchangerInfoManagerEntitiesBase &_base) : base(_base){};

   std::vector<const AOMD::mEntity *> getInfo(int sendto) const
   {
      std::vector<const AOMD::mEntity *> val = std::move(base.created_entities_from[sendto]);
      base.created_entities_from.erase(sendto);
      return val;
   }

   void setInfo(const information_key_t &key, information_t &info, int recv_from)
   {
      AOMD::mEntity &entity = const_cast<AOMD::mEntity &>(*key);
      auto po = base.part_man.getPartitionObject(entity);
      po.insert(recv_from, info);
   }

  private:
   exchangerInfoManagerEntitiesBase &base;
};

class exchangerKeyManagerSendBackNewRemoteObjectWhereNeeded
    : public exchangerKeyManagerSimpleSendOrReceiveKeysFromPartitionManagerAOMD
{
  public:
   exchangerKeyManagerSendBackNewRemoteObjectWhereNeeded(const ranksDataManager &_rank_manager, const partitionManager &part_man)
       : rank_manager(_rank_manager)
   {
      MPI_Comm_rank(part_man.getComm(), &rank);
   }

   std::set<int> getMessageRanks(information_key_t lk) const override
   {
      std::set<int> targets;
      const AOMD::mEntity &entity = *lk;
      const ranksData *prank_data = rank_manager.getData(entity);
      if (prank_data)
      {
         const std::vector<int> &old_rank = prank_data->getOldRanks();
         if (!old_rank.empty())
         {
            const int oldowner = old_rank.front();
            if (oldowner != rank) throw;
            const std::vector<int> new_rank = prank_data->getNewRanks();
            const int new_rank_size = new_rank.size();
            if (new_rank_size)
            {
               for (int i : old_rank)
               {
                  if (i != rank) targets.insert(i);
               }
            }
            if (new_rank_size > 1)
            {
               for (int i : new_rank)
               {
                  targets.insert(i);
               }
            }
         }
      }
      return targets;
   }

  private:
   const ranksDataManager &rank_manager;
   int rank;
};

class exchangerInfoManagerSendBackNewRemoteObjectWhereNeeded
{
  public:
   typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
   typedef xtool::send_only_keys_communication_trait communication_trait;
   typedef const AOMD::mEntity *information_key_t;

   exchangerInfoManagerSendBackNewRemoteObjectWhereNeeded(exchangerInfoManagerEntitiesBase &_base) : base(_base){};

   size_t getApproxDataSize() const { return (sizeof(void *) + sizeof(int)) * 3; };

   information_key_t localObjectKey(const AOMD::mEntity &lo) const { return &lo; };

   void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto) const
   {
      const AOMD::mEntity &entity = *key;
      const ranksData *prank_data = base.rank_manager.getData(entity);
      assert(prank_data);
      const auto po = base.part_man.getConstPartitionObject(entity);
      AOMD::mEntity *pentity_on_sendto = const_cast<AOMD::mEntity *>(po.getRemoteObjectOn(sendto));
      assert(pentity_on_sendto);
      std::vector<int> new_rank = prank_data->getNewRanks();
      auto itend = std::remove_if(new_rank.begin(), new_rank.end(), [&sendto](const int &i) { return (i == sendto); });
      new_rank.erase(itend, new_rank.end());
      const int nbro_to_send = new_rank.size();
      assert(nbro_to_send > 0);
      buff.pack(&pentity_on_sendto, 1, MPI_AINT);
      buff.pack(&nbro_to_send, 1, MPI_INT);
      std::vector<exchangerInfoManagerEntitiesBase::ro_t> ro_buff;
      ro_buff.reserve(nbro_to_send);
      for (int rid : new_rank)
      {
         AOMD::mEntity *prentity_address = const_cast<AOMD::mEntity *>(po.getRemoteObjectOn(rid));
         assert(prentity_address);
         ro_buff.emplace_back(exchangerInfoManagerEntitiesBase::ro_t{rid, prentity_address});
      }
      buff.pack(ro_buff.data(), nbro_to_send, base.mpi_ro_t);
   }

   void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom)
   {
      while (!buff.exhausted())
      {
         AOMD::mEntity *pentity = nullptr;
         buff.unPack(&pentity, 1, MPI_AINT);
         assert(pentity);
         AOMD::mEntity &entity = *pentity;
         base.unPackRemotesEntities(buff, entity);
      }
   }

  private:
   exchangerInfoManagerEntitiesBase &base;
};

class exchangerKeyManagerDeleteEntities : public exchangerKeyManagerSimpleSendOrReceiveKeysFromPartitionManagerAOMD
{
  public:
   exchangerKeyManagerDeleteEntities(const partitionManager &_part_man, const ranksDataManager &_rank_manager)
       : part_man(_part_man), rank_manager(_rank_manager)
   {
   }

   std::set<int> getMessageRanks(information_key_t lk) const override
   {
      std::set<int> targets;
      int rank;
      MPI_Comm_rank(part_man.getComm(), &rank);
      const AOMD::mEntity &entity = *lk;
      const auto *prank_data_entity = rank_manager.getData(entity);
      if (prank_data_entity)
      {
         // I'am the old owner
         // get procs where  copy of entity need to be deleted.
         std::vector<int> ranks_to_remove = prank_data_entity->getRemoveRanks();
         if (!ranks_to_remove.empty())
         {
            targets.insert(rank);
            auto cpo_entity = part_man.getConstPartitionObject(entity);
            for (auto ro_entity : cpo_entity.getRemoteObjectsCollectionRange())
            {
               targets.insert(ro_entity.getProcId());
            }
         }
      }
      return targets;
   }

  private:
   const partitionManager &part_man;
   const ranksDataManager &rank_manager;
};

class exchangerInfoManagerDeleteEntities
{
  public:
   typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
   typedef xtool::send_only_keys_communication_trait communication_trait;
   typedef const AOMD::mEntity *information_key_t;

   exchangerInfoManagerDeleteEntities(exchangerInfoManagerEntitiesBase &_base) : base(_base){};

   size_t getApproxDataSize() const { return sizeof(void *) + sizeof(int) * 2; };

   void getInfo(const information_key_t key, xtool::xMpiInputBuffer &buff, int sendto) const
   {
      const AOMD::mEntity &entity = *key;
      const auto *prank_data_entity = base.rank_manager.getData(entity);
      assert(prank_data_entity);
      std::vector<int> ranks_to_remove = prank_data_entity->getRemoveRanks();
      assert(!ranks_to_remove.empty());
      auto po_entity = base.part_man.getConstPartitionObject(entity);
      const AOMD::mEntity *pr_entity = po_entity.getRemoteObjectOn(sendto);
      buff.pack(&pr_entity, 1, MPI_AINT);
      const int ranks_to_remove_size = ranks_to_remove.size();
      buff.pack(&ranks_to_remove_size, 1, MPI_INT);
      buff.pack(ranks_to_remove.data(), ranks_to_remove_size, MPI_INT);
   }

   void setInfo(const xtool::xMpiOutputBuffer &buff, int receivedfrom)
   {
      int rank;
      MPI_Comm_rank(base.part_man.getComm(), &rank);
      // std::cout << "Set  info, send to " << receivedfrom << std::endl;
      while (!buff.exhausted())
      {
         AOMD::mEntity *pentity;
         buff.unPack(&pentity, 1, MPI_AINT);
         int ranks_to_remove_size;
         buff.unPack(&ranks_to_remove_size, 1, MPI_INT);
         std::vector<int> ranks_to_remove(ranks_to_remove_size);
         buff.unPack(ranks_to_remove.data(), ranks_to_remove_size, MPI_INT);
         auto po_entity = base.part_man.getPartitionObject(*pentity);
         for (int rrank : ranks_to_remove)
         {
            assert(po_entity.getRemoteObjectOn(rrank));
            po_entity.remove(rrank);
         }
         auto itrank = std::find_if(ranks_to_remove.begin(), ranks_to_remove.end(), [&rank](int rrank) { return rrank == rank; });
         if (itrank != ranks_to_remove.end())
         {
            base.delMoreDataEntity(*pentity);
            base.part_man.remove(*pentity);
            base.rank_manager.deleteData(*pentity);
            base.m.DEL_updateAdj(pentity);
         }
      }
   }

  private:
   exchangerInfoManagerEntitiesBase &base;
};

void migrateEntities(AOMD::mMesh &m, partitionManager &part_man, xinterface::aomd::xAttachedDataManagerAOMD<int> &targetproc,
                     moreDataEntity *pmore_data_entity, bool fixed_vertex_id)
{
   typedef const AOMD::mEntity *key_t;
   typedef xtool::xKeyContainerSendOrRecv<key_t> keyContainer;
   // set the rank_manager
   ranksDataManager rank_manager = setRankManager(m, part_man, targetproc);
   targetproc.clear();

   keyContainer exchanger_key_container(part_man.getComm());

   exchangerKeyManagerCreateEntities exchanger_key_manager_create_entities(part_man, rank_manager);
   exchangerKeyManagerSendBackNewRemoteObjectWhereNeeded exchanger_key_manager_sbnro(rank_manager, part_man);
   exchangerKeyManagerDeleteEntities exchanger_key_manager_de(part_man, rank_manager);

   exchangerInfoManagerEntitiesBase exchanger_info_manager_base(m, part_man, rank_manager, pmore_data_entity);
   exchangerInfoManagerCreateVertices exchanger_info_manager_create_vertices(exchanger_info_manager_base, fixed_vertex_id);
   exchangerInfoManagerCreateEdges exchanger_info_manager_create_edges(exchanger_info_manager_base);
   exchangerInfoManagerCreateTriangles exchanger_info_manager_create_triangles(exchanger_info_manager_base);
   exchangerInfoManagerCreateTetrahedra exchanger_info_manager_create_tetrahedra(exchanger_info_manager_base);
   std::array<exchangerInfoManagerCreateEntities *, 4> exchanger_info_manager_create_level = {
       {&exchanger_info_manager_create_vertices, &exchanger_info_manager_create_edges, &exchanger_info_manager_create_triangles,
        &exchanger_info_manager_create_tetrahedra}};
   exchangerInfoManagerSendBackEntitiesAddresses exchanger_info_manager_send_back_entities_addresses(exchanger_info_manager_base);
   exchangerInfoManagerSendBackNewRemoteObjectWhereNeeded exchanger_info_manager_sbnro(exchanger_info_manager_base);
   exchangerInfoManagerDeleteEntities exchanger_info_manager_de(exchanger_info_manager_base);
   // Creating new entities and update partition manager
   for (size_t dim = 0; dim <= 3; ++dim)
   {
      exchanger_key_container.accumulateKeys(m.begin(dim), m.end(dim), exchanger_key_manager_create_entities);
      exchangeInformation(exchanger_key_container, *exchanger_info_manager_create_level[dim]);
      exchangeInformation(exchanger_key_container, exchanger_info_manager_send_back_entities_addresses);
      exchanger_key_container.clearKeys();
      exchanger_key_container.accumulateKeys(m.begin(dim), m.end(dim), exchanger_key_manager_sbnro);
      exchangeInformation(exchanger_key_container, exchanger_info_manager_sbnro);
      exchanger_key_container.clearKeys();
   }
   // Cleaning up
   for (int dim = 3; dim >= 0; --dim)
   {
      exchanger_key_container.accumulateKeys(m.begin(dim), m.end(dim), exchanger_key_manager_de);
      exchangeInformation(exchanger_key_container, exchanger_info_manager_de);
      exchanger_key_container.clearKeys();
   }
};

xinterface::aomd::xAttachedDataManagerAOMD<int> setRandomPartition(AOMD::mMesh &m, int nb_part)
{
   assert(nb_part > 0);
   std::function<int()> part;
   if (nb_part == 1)
      part = []() { return 0; };
   else
      part = [nb_part]() { return rand() % nb_part; };
   xinterface::aomd::xAttachedDataManagerAOMD<int> targetproc;
   for (auto pr : range(m, 3)) targetproc.setData(*pr) = part();
   for (int dim = 2; dim >= 0; --dim)
   {
      for (auto pent : range(m, dim))
      {
         if (!pent->size(dim + 1)) targetproc.setData(*pent) = part();
      }
   }
   return targetproc;
}

xinterface::aomd::xAttachedDataManagerAOMD<int> setParmetisPartition(AOMD::mMesh &m, const partitionManager &part_man,
                                                                     int nb_part)
{
   xinterface::aomd::xAttachedDataManagerAOMD<int> global_vertex_id = globalVerticesNumbering(m, part_man);
   return setParmetisPartition(m, part_man, global_vertex_id, nb_part);
}

xinterface::aomd::xAttachedDataManagerAOMD<int> setParmetisPartition(
    AOMD::mMesh &m, const partitionManager &part_man, const xinterface::aomd::xAttachedDataManagerAOMD<int> &global_vertex_id,
    int nb_part)
{
   const auto rrange = range(m, 3);
   const auto frange = range(m, 2);
   const auto erange = range(m, 1);
   const auto vrange = range(m, 0);

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
      eptr[0] = 0;
      int pos = 1;
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
            }
         }
      }

      xinterface::parmetis::ParMetisInterface::parmetis_indx_t wgtflag = 0;
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t numflag = 0;
      xinterface::parmetis::ParMetisInterface::parmetis_indx_t ncon = 1;
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
      xinterface::parmetis::ParMetisInterface::PartMeshKway(elmdist.data(), eptr.data(), eind.data(), nullptr, wgtflag, numflag,
                                                            ncon, ncommonnodes, nparts, tpwgts.data(), ubvec.data(),
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

void exportGMSH_V2_Dist(const AOMD::mMesh &m, const partitionManager &part_man, const std::string &filename_base)
{
   int mpi_size, mpi_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

   std::string outfilename = filename_base + "_NP" + std::to_string(mpi_size) + "_P" + std::to_string(mpi_rank) + ".msh";
   // std::cout << "EXPORTing " << outfilename << std::endl;
   std::ofstream out(outfilename.c_str());
   out << "$MeshFormat" << std::endl;
   out << "2.2 0 " << sizeof(double) << std::endl;
   out << "$EndMeshFormat" << std::endl;
   out << "$Nodes" << std::endl;
   out << m.size(0) << std::endl;
   int nbelem = 0;
   for (auto pv : range(m, 0))
   {
      const AOMD::mVertex &v = static_cast<AOMD::mVertex &>(*pv);
      const Trellis_Util::mPoint p = v.point();
      out << v.getId() << " " << p(0) << " " << p(1) << " " << p(2) << "\n";
      if (v.getClassification()->dim() == 0)
         ++nbelem;
      else if (part_man.hasRemoteCopy(*pv))
         ++nbelem;
   }
   out << "$EndNodes" << std::endl;

   for (int dim = 1; dim <= 2; ++dim)
   {
      for (auto pe : range(m, dim))
      {
         const AOMD::mEntity &e = static_cast<AOMD::mEntity &>(*pe);
         if (e.getClassification()->dim() == dim)
            ++nbelem;
         else if (part_man.hasRemoteCopy(*pe))
            ++nbelem;
      }
   }

   nbelem += m.size(3);
   out << "$Elements" << std::endl;
   out << nbelem << std::endl;
   int k = 1;
   for (auto pv : range(m, 0))
   {
      AOMD::mVertex &vertex = static_cast<AOMD::mVertex &>(*pv);
      if ((vertex.getClassification()->dim() == 0) || part_man.hasRemoteCopy(vertex))
      {
         const int tag = vertex.getClassification()->tag();
         const int id0 = vertex.getId();
         out << k++ << " 15 4 " << tag << " " << tag << " 1 " << mpi_rank + 1 << " " << id0 << "\n";
      }
   }
   for (auto pe : range(m, 1))
   {
      AOMD::mEntity &edge = *pe;
      if ((edge.getClassification()->dim() == 1) || part_man.hasRemoteCopy(edge))
      {
         const int tag = edge.getClassification()->tag();
         const int id0 = edge.get(0, 0)->getId();
         const int id1 = edge.get(0, 1)->getId();
         out << k++ << " 1 4 " << tag << " " << tag << " 1 " << mpi_rank + 1 << " " << id0 << " " << id1 << "\n";
      }
   }
   for (auto pf : range(m, 2))
   {
      AOMD::mEntity &tri = *pf;
      if ((tri.getClassification()->dim() == 2) || part_man.hasRemoteCopy(tri))
      {
         const int tag = tri.getClassification()->tag();
         const int id0 = tri.get(0, 0)->getId();
         const int id1 = tri.get(0, 1)->getId();
         const int id2 = tri.get(0, 2)->getId();

         out << k++ << " 2 4 " << tag << " " << tag << " 1 " << mpi_rank + 1 << " " << id0 << " " << id1 << " " << id2 << "\n";
      }
   }
   for (auto pr : range(m, 3))
   {
      AOMD::mEntity &tet = *pr;
      const int tag = tet.getClassification()->tag();
      const int id0 = tet.get(0, 0)->getId();
      const int id1 = tet.get(0, 1)->getId();
      const int id2 = tet.get(0, 2)->getId();
      const int id3 = tet.get(0, 3)->getId();

      out << k++ << " 4 4 " << tag << " " << tag << " 1 " << mpi_rank + 1 << " " << id0 << " " << id1 << " " << id2 << " " << id3
          << "\n";
   }
   out << "$EndElements" << std::endl;
   out.close();
}

void exportPartition(AOMD::mMesh &m, std::ofstream &f,
                     xinterface::aomd::xAttachedDataManagerAOMD<int> &targetprocfortoplevelentities,
                     xinterface::aomd::xAttachedDataManagerAOMD<int> &entities_id)
{
   for (auto e : xtool::make_range(targetprocfortoplevelentities.beginKey(), targetprocfortoplevelentities.endKey()))
   {
      int *p = entities_id.getData(*e);
      if (p)
      {
         f << *p << " " << *(targetprocfortoplevelentities.getData(*e)) << endl;
      }
      else
      {
         cout << "We do not know entity id with AOMD id " << e->getId() << endl;
         throw -234;
      }
   }
}

}  // namespace xmeshtool
