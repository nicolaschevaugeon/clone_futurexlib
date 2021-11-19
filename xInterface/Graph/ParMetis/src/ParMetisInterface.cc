/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/
// std
#include <iostream>
#include <vector>
#include <cassert>

#include "ParMetisInterface.h"
#include "parmetis.h"

namespace xinterface
{
  namespace parmetis
  {


 
    void ParMetisInterface::PartKway(const parmetis_indx_t *vtxdist, const parmetis_indx_t *xadj, const parmetis_indx_t *adjncy,
				     const parmetis_indx_t * vwgt, const parmetis_indx_t * adjwgt, const parmetis_indx_t &wgtflag,
				     const parmetis_indx_t & numflag, const parmetis_indx_t & ncon, const parmetis_indx_t & nparts,
				     const parmetis_real_t *tpwgts, const parmetis_real_t *ubvec,
				     const parmetis_indx_t *options, MPI_Comm world,
				     parmetis_indx_t &edgecut,parmetis_indx_t *part)
    {
      int res=0;
#if ( PARMETIS_MAJOR_VERSION > 3 )
      res =
#endif
	ParMETIS_V3_PartKway (
			      const_cast < parmetis_indx_t * >( vtxdist )
			      , const_cast < parmetis_indx_t * >( xadj )
			      , const_cast < parmetis_indx_t * >( adjncy )
			      , const_cast < parmetis_indx_t * >( vwgt )
			      , const_cast < parmetis_indx_t * >( adjwgt )
			      , const_cast < parmetis_indx_t * >( &wgtflag )
			      , const_cast < parmetis_indx_t * >( &numflag )
			      , const_cast < parmetis_indx_t * >( &ncon )
			      , const_cast < parmetis_indx_t * >( &nparts )
			      , const_cast < parmetis_real_t * >( tpwgts )
			      , const_cast < parmetis_real_t * >( ubvec )
			      , const_cast < parmetis_indx_t * >( options )
			      , &edgecut
			      , part
			      , const_cast < MPI_Comm * >( &world ));
      error_treatement(res);
      return;
    }

    void ParMetisInterface::PartKway(const parmetis_indx_t *vtxdist, const parmetis_indx_t *xadj, const parmetis_indx_t *adjncy,
				     const parmetis_indx_t * vwgt, const parmetis_indx_t * adjwgt, const parmetis_indx_t &wgtflag,
				     const parmetis_indx_t & ncon, const parmetis_indx_t & nparts, MPI_Comm world,
				     parmetis_indx_t &edgecut,parmetis_indx_t *part)
    {
      std::vector <  parmetis_real_t > tpwgts(nparts*ncon,1.f/nparts);
      std::vector <  parmetis_real_t > ubvec(ncon,1.05f);
      parmetis_indx_t numflag = 0;
      parmetis_indx_t options[3] = {0,0,15};
      int res=0;
#if ( PARMETIS_MAJOR_VERSION > 3 )
      res =
#endif
	ParMETIS_V3_PartKway (
			      const_cast < parmetis_indx_t * >( vtxdist )
			      , const_cast < parmetis_indx_t * >( xadj )
			      , const_cast < parmetis_indx_t * >( adjncy )
			      , const_cast < parmetis_indx_t * >( vwgt )
			      , const_cast < parmetis_indx_t * >( adjwgt )
			      , const_cast < parmetis_indx_t * >( &wgtflag )
			      , &numflag
			      , const_cast < parmetis_indx_t * >( &ncon )
			      , const_cast < parmetis_indx_t * >( &nparts )
			      , &tpwgts[0]
			      , &ubvec[0]
			      , &options[0]
			      , &edgecut
			      , part
			      , const_cast < MPI_Comm * >( &world ));
      error_treatement(res);
      return;

    }
    void ParMetisInterface::PartMeshKway(const parmetis_indx_t *elmdist, const parmetis_indx_t *eptr, const parmetis_indx_t *eind,
					 const parmetis_indx_t *elmwgt, const parmetis_indx_t  &wgtflag, const parmetis_indx_t &numflag,
					 const parmetis_indx_t &ncon, const parmetis_indx_t &ncommonnodes, const parmetis_indx_t &nparts,
					 const parmetis_real_t *tpwgts, const parmetis_real_t *ubvec,
					 const parmetis_indx_t *options, MPI_Comm world,
					 parmetis_indx_t &edgecut,parmetis_indx_t *part)
    {
      int res=0;
#if ( PARMETIS_MAJOR_VERSION > 3 )
      res =
#endif
	ParMETIS_V3_PartMeshKway (
				  const_cast < parmetis_indx_t * >  ( elmdist )
				  , const_cast < parmetis_indx_t * >( eptr )
				  , const_cast < parmetis_indx_t * >( eind )
				  , const_cast < parmetis_indx_t * >( elmwgt )
				  , const_cast < parmetis_indx_t * >( &wgtflag )
				  , const_cast < parmetis_indx_t * >( &numflag )
				  , const_cast < parmetis_indx_t * >( &ncon )
				  , const_cast < parmetis_indx_t * >( &ncommonnodes )
				  , const_cast < parmetis_indx_t * >( &nparts )
				  , const_cast < parmetis_real_t * >( tpwgts )
				  , const_cast < parmetis_real_t * >( ubvec )
				  , const_cast < parmetis_indx_t * >( options )
				  , &edgecut
				  , part
				  , const_cast < MPI_Comm * >( &world ));
      error_treatement(res);
      return;
    }
    void ParMetisInterface::error_treatement(int res)
    {
#if ( PARMETIS_MAJOR_VERSION > 3 )
      if (res != METIS_OK)
	{
	  std::cout << "Error in calling ParMetis function" << std::endl;
	  switch (res)
	    {
            case METIS_ERROR_INPUT :
	      {
		std::cout << "Returned due to erroneous inputs and/or options" << std::endl;
		break;
	      }
            case  METIS_ERROR_MEMORY :
	      {
		std::cout << "Returned due to insufficient memory " << std::endl;
		break;
	      }
            default :
            case  METIS_ERROR :
	      {
		std::cout << "Some unducomented error" << std::endl;
		break;
	      }
	    }
	  throw res;
	}
#endif
    }

  } // end subnamespace

} // end namespace
