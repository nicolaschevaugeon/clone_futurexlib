/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/
#ifndef XSPLITAOMD_IMP_H
#define XSPLITAOMD_IMP_H


namespace xinterface
{

  namespace xtemplaterefinemesh 
  {

    template < typename CRIT >
      xMeshSplitUserCriteriaAOMD < CRIT >::xMeshSplitUserCriteriaAOMD(CRIT &criteria_, int max_it_split_ ) :
    splitor(criteria_,max_it_split_)
      {
	registerTransfer(&trans_classification);
      }
    template < typename CRIT >
      inline void xMeshSplitUserCriteriaAOMD < CRIT >::registerTransfer(xSplitTransferBaseAOMD* trans_information)
      {
	splitor.registerTransfer(trans_information);
      }
    template < typename CRIT >
      inline bool xMeshSplitUserCriteriaAOMD < CRIT >::split(AOMD::mMesh &mesh,int target_dim_, partman_t & part_man_)
      {
	xMeshSplitWarperAOMD warper(mesh);
	return splitor.split(warper,target_dim_,part_man_);
      }
    template < typename CRIT >
      inline void xMeshSplitUserCriteriaAOMD < CRIT >::clearInternal()
      {
	splitor.clearInternal();
      }

  } // end subnamespace

} // end namespace

#endif
