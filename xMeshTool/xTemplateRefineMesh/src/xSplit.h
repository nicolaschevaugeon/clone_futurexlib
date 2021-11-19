/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE & LICENSE files for terms 
   and conditions.
*/

#ifndef XSPLIT_H
#define XSPLIT_H

// std
#include <iostream>
#include <functional>
#include <string>
#include <set>
#include <array>
#include <list>
#include <vector>
#include <forward_list>


// eXlibris_tools
#include "xPartitionManager.h"
#include "xDataExchanger.h"

// xmeshtool
#include "xSplitInfoManager.h"


// parallel
#include "mpi.h"


namespace xmeshtool
{


  template < typename T, template < typename > class DM >
    class xMeshSplitUserCriteria
    {
    public:
      /// type from traits T
      typedef typename T::criteria_type criteria_t;
      typedef typename T::transfert_info_type transfert_info_t;
      typedef typename T::warper_type warper_t;
      typedef typename warper_t::entity_type entity_t;



      /// the type of partition manager used by xMeshSplitUserCriteria
      typedef xtool::xPartitionManager < DM > partman_t;

      /// Constructor
      //! \param criteria_  : criteria that tels if a entity have to be split or not.
      //!                      It's user responsibility to provide an object that have at least
      //!                      a () operator which take a T reference as argument and return a boolean to
      //!                      say if this entity is to be split or not.
      //!                      It must also have a update method wich based on new
      //!                      state of the mesh, update eventually some information to give a correct answer fo () operator
      //
      //!                      In future we hope to have an extra method to unsplit ...
      //
      xMeshSplitUserCriteria(criteria_t &criteria_, int max_it_split_ = 30);
      /// registration
      //! \param  trans_information : this is a transferring mechanism. It provide at least :
      //!                                  void collect(AOMD::mEntity *e) : a method which collect informations related to e
      //!                                  void transfer(AOMD::mEntity *se) : a method which transfer informations previously collected
      //!                                    and set them in se.
      //!                                 It is used when a entity (e) is split in new entities (se)
      //
      void registerTransfer(transfert_info_t* trans_information);
      bool split(warper_t & warp,int target_dim_, partman_t & part_man_);
      void clearInternal();
    private:
      // private data
      typedef xtool::xRange < typename warper_t::iter_type > range_t;
      typedef std::array < entity_t *, 3 >  t1_t;
      typedef std::array < entity_t *, 5 >  t2_t;
      std::set < entity_t * > edge_face_to_be_treated;
      warper_t * warper;
      criteria_t & criteria;
      std::forward_list < transfert_info_t * > trans_information_container;
      int target_dim,max_it_split;
      DM < entity_t * > hanging_down,hanging_up;
      DM < t1_t > hanging_data_t1;
      DM < t2_t > hanging_data_t2;
      partman_t *part_man;
      MPI_Comm *world;

      // tet conectivity chosed for this class (Note : comes from AOMD)
      const int Tev[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
      const int Tfv[4][3] = {{0,1,2},{0,1,3},{1,2,3},{0,2,3}};
      const int Tfe[4][3] = {{0,1,2},{0,4,3},{1,5,4},{2,3,5}};
      const int Tve[4][4] = {{-1,0,2,3},{0,-1,1,4},{2,1,-1,5},{3,4,5,-1}};

      // All key/info manager/container are suposed to use the same information key type. So we take one arbitrary
      typedef  typename splitKeyManagerSendOnly < partman_t, entity_t >::information_key_t information_key_t;

      //  splitInfoManagerNewHangingTreatement : frontier element says to its remote that it have some new hanging frontier element
      xtool::xKeyContainerSendOrRecv < information_key_t > *new_hanging_key_cont;
      // key manager for Send Only pattern
      splitKeyManagerSendOnly<partman_t, entity_t > *kmso;
      // key manager for Send And receive pattern
      splitKeyManagerSendAndReceive<partman_t, entity_t > *kmsar;

      // private function object to private methode
      std::function < void(entity_t *) > select_hanging_edge;
      std::function < void(entity_t *) > select_hanging_face;
      std::function < bool(entity_t &,entity_t &) > delete_hanging_face_info;
      std::function < bool(entity_t &,entity_t &) > delete_hanging_edge_info_from_tet;
      std::function < bool(entity_t &,entity_t &) > delete_hanging_edge_info_from_face;
      std::function < void(entity_t &) > treat_new_hanging;
      std::function < entity_t *(entity_t &) > split_elem;
      std::function < entity_t *(entity_t &) > split_related;

      // private method
      // ====================================
      void checkOneLevelCriteria( entity_t *e, std::set < entity_t * > &to_split_onelevel, std::set < entity_t * > &bnd_forced_to_split);
      void addOneLevelCriteria(std::set < entity_t * > &bnd_which_add_to_to_split, std::list < entity_t * > &to_split);

      // spliter
      entity_t * splitEdge(entity_t &edge);
      entity_t * splitFace(entity_t &face);
      entity_t * splitTet(entity_t &tet);
      entity_t * splitTetRelatedToHanging(entity_t &tet);

      // hanging cleanner
      bool deleteFaceWithCondition(entity_t &,entity_t &);
      bool deleteEdgeWithCondition(entity_t &,entity_t &);
      bool dontDeleteEntityWithCondition(entity_t &,entity_t &);
      void deleteHangingInfoForFace(entity_t &,entity_t &);
      void deleteHangingInfoForEdge(entity_t &,entity_t &);

      // create
      entity_t * createTetFromVertex(entity_t *,entity_t *,entity_t *,entity_t *, bool *);


      //hanging intermediate synchronisation
      void hangingSynchro(entity_t &);
      void removeDangling();

      // utility
      void dontSelect(entity_t *);
      void doSelect(entity_t *);
      void collect(entity_t *e);
      void transfer(entity_t *e);
      const entity_t *smallestAddress(entity_t *e);

    };



} // end namespace

#include "xSplit_imp.h"

#endif
