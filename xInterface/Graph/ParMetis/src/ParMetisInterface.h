/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/
#ifndef PARMETISINTERFACE_H
#define PARMETISINTERFACE_H
#include "mpi.h"

 
namespace xinterface
{
  namespace parmetis
  {
    class DistCSRGraphParMetisbase
    {};
    class ParMetisParameterBase
    {};

    class ParMetisInterface
    {

    public:
      // types
#include "ParMetisInterface_types.h"

      // In most case those function return 1 if ok (last parmetis version). For older version ok may be one but function
      // did not work properly.

      //! partition a graph giving its vertex and edges
      //! most general
      static void PartKway( const parmetis_indx_t *vtxdist, const parmetis_indx_t *xadj, const parmetis_indx_t *adjncy,
			    const parmetis_indx_t * vwgt, const parmetis_indx_t * adjwgt, const parmetis_indx_t &wgtflag,
			    const parmetis_indx_t & numflag, const parmetis_indx_t & ncon, const parmetis_indx_t & nparts,
			    const parmetis_real_t *tpwgts, const parmetis_real_t *ubvec,
			    const parmetis_indx_t *options, MPI_Comm world,
			    parmetis_indx_t &edgecut,parmetis_indx_t *part);
      /*! partition a graph giving its vertex and edges
       * less general:
       *  - numflag defaulted to 0
       *  - tpwgts defaulted to 1/nb_proc
       *  - ubvec defaulted to 1.05
       *  - option defaulted to not extra output
       */
      static void  PartKway(const parmetis_indx_t *vtxdist, const parmetis_indx_t *xadj, const parmetis_indx_t *adjncy,
			    const parmetis_indx_t * vwgt, const parmetis_indx_t * adjwgt, const parmetis_indx_t &wgtflag,
			    const parmetis_indx_t & ncon, const parmetis_indx_t & nparts, MPI_Comm world,
			    parmetis_indx_t &edgecut,parmetis_indx_t *part);

      // partition a graph giving its vertex and edges
      // with graph and parameter object TO BE DONE if ever
      void PartKway(const DistCSRGraphParMetisbase &graph, const ParMetisParameterBase *param){; /* to be done */}


      /*! partition a mesh generating its dual graph to compute the partition
       */
      static void PartMeshKway(const parmetis_indx_t *elmdist, const parmetis_indx_t *eptr, const parmetis_indx_t *eind,
			       const parmetis_indx_t *elmwgt, const parmetis_indx_t  &wgtflag, const parmetis_indx_t &numflag,
			       const parmetis_indx_t &ncon, const parmetis_indx_t &ncommonnodes, const parmetis_indx_t &nparts,
			       const parmetis_real_t *tpwgts, const parmetis_real_t *ubvec,
			       const parmetis_indx_t *options, MPI_Comm world,
			       parmetis_indx_t &edgecut,parmetis_indx_t *part);

    private:
      static void error_treatement(int res);

    };

    // to be done
    class DistCSRGraphParMetis : public DistCSRGraphParMetisbase
    {
    public:
      //void connect(ParMetisInterface::parmetis_indx_t **vtxdist,....
    private:
      ParMetisInterface::parmetis_indx_t *vtxdist;
      ParMetisInterface::parmetis_indx_t *xadj;
      ParMetisInterface::parmetis_indx_t *adjncy;
      ParMetisInterface::parmetis_indx_t * vwgt;
      ParMetisInterface::parmetis_indx_t * adjwgt;
      MPI_Comm world;
    };
    class ParMetisParameter : public ParMetisParameterBase
    {
      /*
	ParMetisInterface::parmetis_indx_t &wgtflag;

	ParMetisInterface::parmetis_indx_t & numflag;
	ParMetisInterface::parmetis_indx_t & ncon;
	ParMetisInterface::parmetis_indx_t & nparts;
      */

      ParMetisInterface::parmetis_real_t *tpwgts;
      ParMetisInterface::parmetis_real_t *ubvec;

    };



  } // end subnamespace
} // end namespace
#endif
