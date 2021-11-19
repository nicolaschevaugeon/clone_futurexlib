/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/
#ifndef XSPLITAOMD_H
#define XSPLITAOMD_H

// AOMD
#include "mMesh.h"

// xtool
#include "xPartitionManager.h"

// xAOMDInterfaceGeneral
#include "xAttachedDataManagerAOMD.h"

// xmeshtool
#include "xSplit.h"

namespace xinterface
{

  namespace xtemplaterefinemesh
  {
    class xSplitTransferBaseAOMD
    {
    public:
      virtual void collect(AOMD::mEntity *e) = 0;
      virtual void transfer(AOMD::mEntity *se) = 0;
    };
    class xSplitTransferAOMDGeom : public xSplitTransferBaseAOMD
    {
    public:
      ///  A method which collect geometric classification related to e
      void collect(AOMD::mEntity *e) override;
      /// A method which transfer classification previously collected
      void transfer(AOMD::mEntity *se) override;
    private:
      AOMD::GEntity *g;
    };

    // elmentary function warper
    class xMeshSplitWarperAOMD
    {
    public:
      typedef AOMD::mMesh mesh_type;
      typedef AOMD::mEntity entity_type;
      typedef typename AOMD::mMesh::iter iter_type;
      xMeshSplitWarperAOMD(AOMD::mMesh &m_);
      AOMD::mEntity * createCOGVertex        (AOMD::mEntity &edge);
      AOMD::mEntity * createEdge             (AOMD::mEntity *v1, AOMD::mEntity *v2);
      AOMD::mEntity * createFaceWithEdge     (AOMD::mEntity *e1, AOMD::mEntity *e2, AOMD::mEntity *e3);
      AOMD::mEntity * createTetWithFace      (AOMD::mEntity *f1, AOMD::mEntity *f2, AOMD::mEntity *f3, AOMD::mEntity *f4);
      AOMD::mEntity * getEdgeFromVertex      (AOMD::mEntity *v1, AOMD::mEntity *v2);
      AOMD::mEntity * getFaceFromVertex      (AOMD::mEntity *v1, AOMD::mEntity *v2,AOMD::mEntity *v3);
      void            removeEntityWithAdj    (AOMD::mEntity &e);
      void            printEntity            (const AOMD::mEntity *e, bool deep);
      void            printMesh              (bool deep);
      void            writeToFile            (const std::string &file_name);
      iter_type       begin                  (int what);
      iter_type       end                    (int what);

      // Aspect ratio function 
      double          tetAspectRatio         (AOMD::mEntity *tet);
      double          tetAspectRatio         (AOMD::mEntity *v1, AOMD::mEntity *v2,AOMD::mEntity *v3, AOMD::mEntity *v4);

      // acumulate points
      void            accumulatePointCoord   (const AOMD::mEntity *e);
      void            resetPointCoord        ();
      AOMD::mEntity * createMeanPointCoord   ();

    private:
      AOMD::mMesh &mesh;
      double pc[3];
      int nb_pc;

    };

    template < typename CRIT >
      class xMeshSplitUserCriteriaAOMD
      {

      public:
        typedef xtool::xPartitionManager < xinterface::aomd::xAttachedDataManagerAOMD > partman_t;

        /// Constructor
        //! \param criteria_  : criteria that tels if a entity have to be split or not.
        //!                      It's user responsibility to provide an object that have at least
        //!                      a () operator which take a T reference as argument and return a boolean to
        //                       say if this entity is to be split or not.
        //                       It must also have a update method which take a mesh_type reference and based on new
        //                       state of the mesh update eventually some information to give a correct answer fo () operator
        xMeshSplitUserCriteriaAOMD(CRIT &criteria_, int max_it_split_ = 30);
        /// registration
        //! \param  trans_information : this is a transferring mechanism. It provide at least :
        //!                                  void collect(AOMD::mEntity *e) : a method which collect informations related to e
        //!                                  void transfer(AOMD::mEntity *se) : a method which transfer informations previously collected
        //!                                    and set them in se.
        //!                                 It is used when a entity (e) is split in new entities (se)
        //
        inline void registerTransfer(xSplitTransferBaseAOMD* trans_information);
        bool split(AOMD::mMesh &mesh,int target_dim_, partman_t & part_man_);
        inline void clearInternal();

      private:
        struct SplitAOMDTraits
        {
	  typedef CRIT criteria_type;
	  typedef xSplitTransferBaseAOMD transfert_info_type;
	  typedef xMeshSplitWarperAOMD warper_type;
        };
	xmeshtool::xMeshSplitUserCriteria < SplitAOMDTraits,xinterface::aomd::xAttachedDataManagerAOMD > splitor;
        xSplitTransferAOMDGeom trans_classification;
      };


  } // end subnamespace

} // end namespace

#include "xSplitAOMD_imp.h"

#endif
