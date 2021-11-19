/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef XLOADBALANCETOOLS
#define XLOADBALANCETOOLS
#include <algorithm>
#include <iostream>
#include <vector>

#include "xAttachedDataManagerAOMD.h"
#include "xDataExchanger.h"
#include "xPartitionManager.h"
// graph_tools
#include "ParMetisInterface.h"

namespace xmeshtool
{
class moreDataEntity;
typedef xtool::xPartitionManager<xinterface::aomd::xAttachedDataManagerAOMD> partitionManager;

/// \brief given a mesh, a partition manager and targetproc migrate mesh entity of the mesh accross processor. The call is
/// collective among all the processors.
///
/// The mesh and the partition manager must describe a "valid"  partition : \n
///   - A top level entity (an entity that has no upward adjacencies) must be present only on one proc. \n
///   - The partition manager on entry must be correct. \n
///      Note that if the mesh is empty on all procs but one, the partition manager is correct if empty \n
///      ... It means that you can call this function on a mesh not yet partitionned ...\n
/// The input mesh must : \n
///   - be a mesh of simplices (only TET, TRI, EDGE and VERTEX entities) \n
///   - A complet mesh, meaning it contains all entities:   \n
///      for each tet, all bounding  tri must be present, for each tri all the bounding edges must be present, \n
///      for each edge all vertex must be present   \n
///      Note that it is typically the case in xFem library (see xMesh construction : -> ModifyAllState) \n
///   -  Only the  one level adjacencies are needed : TET<->TRI<->EDGE<->VERTEX \n
///      If more are present, they will be ignored.\n
/// The \param targetproc drive the entities migration. All top level entity of the input mesh must have a value associated to it
/// in this datamanager \n
///     this value tell on which proc one want to have the entity after the call to this function.
/// On output, If successful, on each proc :
///   - the mesh contain only the toplevelentities as defined by target, is complete and has all the one level adjacencies. \n
///     Note that other adjacencies are NOT maintained. \m
///   - The part_man is updated accordingly and correct on each proc. \n
///   pmore_data_entity is an optional argument to let the user communicate and clean up additional data associated to elements
///   data
void migrateEntities(AOMD::mMesh &m, partitionManager &part_man,
                     xinterface::aomd::xAttachedDataManagerAOMD<int> &targetprocfortoplevelentities,
                     moreDataEntity *pmore_data_entity = nullptr, bool fixed_vertex_id = true);

/// \brief moreDataEntity is a base class to define call back used when an entity migrate entities.
///
/// Typical usage is to ensure some data associated to an entity can be retrieved after the migration, or to add some debugging
/// messaga
class moreDataEntity
{
  public:
   /// getInfo is called each time an entity will be send to proc sendto. the user as access to the buffer that will be send and
   /// can pack data into it.
   virtual void getInfo(const AOMD::mEntity &entity, xtool::xMpiInputBuffer &buff, int sendto) const = 0;
   /// setInfo is called just after an entiy has beencreated on a proc. all downward adjancencies to entity are allready present
   /// in the mesh. Data set by the corresponding getInfo in proc with id receivedfrom
   virtual void setInfo(AOMD::mEntity &entity, const xtool::xMpiOutputBuffer &buff, int receivedfrom) = 0;
   /// delEntity is called each time an entity will be delete. typical usage is to clean up some data associated to the entity,
   /// before it disappear from the mesh.
   virtual void delEntity(AOMD::mEntity &entity) = 0;
};

/// \brief setRandomPartition set up random target values between 0 and nb_part for each top level entity in m.
xinterface::aomd::xAttachedDataManagerAOMD<int> setRandomPartition(AOMD::mMesh &m, int nb_part);

/// \brief setParmetisPartition set up target proc values between 0 and nb_part for each top level entity in m, using parmetis. a
/// global numebring is needed for this version.
xinterface::aomd::xAttachedDataManagerAOMD<int> setParmetisPartition(
    AOMD::mMesh &m, const partitionManager &part_man, const xinterface::aomd::xAttachedDataManagerAOMD<int> &global_vertex_id,
    int nb_part);

/// \brief setParmetisPartition set up target proc values between 0 and nb_part for each top level entity in m, using parmetis.
/// Same as above, but the global numbering is created first internally and then destroyed.
xinterface::aomd::xAttachedDataManagerAOMD<int> setParmetisPartition(AOMD::mMesh &m, const partitionManager &part_man,
                                                                     int nb_part);

/// \brief setParmetisPartitionWithWeight set up target proc values between 0 and nb_part for each top level entity in m, using
/// parmetis and given weight.
template <int NBW = 1>
xinterface::aomd::xAttachedDataManagerAOMD<int> setParmetisPartitionWithWeight(
    AOMD::mMesh &m, const partitionManager &part_man, int nb_part,
    xinterface::aomd::xAttachedDataManagerAOMD<std::array<xinterface::parmetis::ParMetisInterface::parmetis_indx_t, NBW>>
        &weights);

/// \brief number all vertices of the mesh so that vertices present in the partitionManager have the same number accros all
/// process. Collective among all the process of the MPI_Comm part_man was build with.
inline xinterface::aomd::xAttachedDataManagerAOMD<int> globalVerticesNumbering(AOMD::mMesh &m, const partitionManager &part_man);

/// \brief number all entities between b and e so that entities present in the partitionManager have the same number accros all
/// process. Collective among all the process of the MPI_Comm part_man was build with.
template <class ITER>
inline xinterface::aomd::xAttachedDataManagerAOMD<int> globalEntitiesNumbering(const ITER &b, const ITER &e,
                                                                               const partitionManager &partman);

/// \brief export partition to a stream
void exportPartition(AOMD::mMesh &m, std::ofstream &f,
                     xinterface::aomd::xAttachedDataManagerAOMD<int> &targetprocfortoplevelentities,
                     xinterface::aomd::xAttachedDataManagerAOMD<int> &entities_id);

void exportGMSH_V2_Dist(const AOMD::mMesh &m, const partitionManager &part_man, const std::string &filename_base);

}  // namespace xmeshtool
#include "xLoadBalanceTools_imp.h"
#endif
