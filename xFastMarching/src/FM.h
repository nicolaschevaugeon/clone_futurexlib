/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FM_h_
#define _FM_h_

#include "FMUtil.h"
#include "FMInternal.h"
#include "FMUpdater.h"


namespace xfastmarching
{
/// Driver for fast marching problem.
/*!
   Compute in ls result of fast marching starting with values at known vertices (given by knownbeg and knownend) ls values at the nodes in this iterator ranges must be known inside ls. "Velocity " or target norm of the gradient at nodes must be given via Ffunc.
   All queries to meshdatabase are done throught the MESHINTERFACE parameter.
   Parameters gls and vn correspond respectively to ls gradient and transported quantity. They are infact kind of optional :  if not given one of the version below is actually launched.
   gls and vn are following the DATAMANAGER concept depicted in FMEntityStorage.h by class entitystorage
   Datas (either gradiant, of VECTOR type, or transported value, of TRANSPORTED type) are associated with vertex (MESHINTERFACE::vertex type) by DATAMANAGER concept.
   !*/
template < template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE, class ITTK, class ITVK, class ITK, class SCAL, class TRANSPORTED >
void fmeik ( const MESHINTERFACE & mi,
             DATAMANAGER < SCAL > & ls,
             const ITTK &knownbeg, const ITTK &knownend,
             const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
             ITK knownit,
             const std::function < SCAL (const typename MESHINTERFACE::vertex &) > & Ffunc,
             SCAL epsilon_ratio,
             DATAMANAGER  <  VECT > & gls,
             DATAMANAGER  < TRANSPORTED > &vn)
{
    internal::FMupdaterGeneric < DATAMANAGER, MESHINTERFACE, VECT, SCAL, TRANSPORTED > updater(mi, ls, Ffunc, epsilon_ratio, gls, vn);
    internal::fmeik(mi, knownbeg, knownend, knownvolatilbeg,knownvolatilend,knownit, updater);
}

/// This version does the fast marching but only compute the gls (levelset gradient at node) additional quantity.
template < template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE, class ITTK, class ITVK, class ITK, class SCAL >
void fmeik ( const MESHINTERFACE & mi,
             DATAMANAGER < SCAL > & ls,
             const ITTK &knownbeg, const ITTK &knownend,
             const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
             ITK knownit,
             const std::function < SCAL (const typename MESHINTERFACE::vertex &) > & Ffunc,
             SCAL epsilon_ratio,
             DATAMANAGER < VECT > & gls)
{
    internal::FMupdaterGeneric < DATAMANAGER, MESHINTERFACE, VECT, SCAL, int > updater(mi, ls, Ffunc, epsilon_ratio, gls);
    internal::fmeik(mi, knownbeg, knownend,knownvolatilbeg,knownvolatilend, knownit, updater);
}

/// This version does the fast marching but without computing any additional quantities.
template < template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE,  class ITTK, class ITVK,  class ITK, class SCAL >
void fmeik ( const MESHINTERFACE & mi,
             DATAMANAGER < SCAL > & ls,
             const ITTK &knownbeg, const ITTK &knownend, const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
             ITK knownit,
             const std::function < SCAL (const typename MESHINTERFACE::vertex &) > & Ffunc,
             SCAL epsilon_ratio)
{
    internal::FMupdaterGeneric < DATAMANAGER, MESHINTERFACE, VECT, SCAL, int > updater(mi, ls, Ffunc, epsilon_ratio );
    internal::fmeik(mi, knownbeg, knownend, knownvolatilbeg,knownvolatilend, knownit, updater);
}

#ifdef HAVE_XGRAPH
/// This version does the fast marching and compute the flux graphs "from" and "to"
template <  template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE,  class ITTK, class ITVK,  class SCAL >
void fmeik ( const MESHINTERFACE &mi,
             DATAMANAGER < SCAL > &ls,
             const ITTK &knownbeg, const ITTK &knownend,
             const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
             const std::function < SCAL (const typename MESHINTERFACE::vertex &) > &Ffunc,
             SCAL epsilon_ratio,
             DATAMANAGER < xgraph::nodeFrom < typename MESHINTERFACE::vertex,SCAL,VECT,3 > > &flux_from,
             DATAMANAGER < xgraph::nodeTo < typename MESHINTERFACE::vertex,SCAL,VECT > > &flux_to
             )
{
    internal::FMupdaterFlux < DATAMANAGER, MESHINTERFACE, VECT, SCAL > updater(mi, ls, Ffunc, epsilon_ratio, flux_from, flux_to);
    internal::fmeikFlux(mi, knownbeg, knownend, knownvolatilbeg,knownvolatilend, updater);
}
/// This version does the fast marching but and compute the flux graph "to"
template <  template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE,  class ITTK, class ITVK,  class SCAL >
void fmeik ( const MESHINTERFACE &mi,
             DATAMANAGER < SCAL > &ls,
             const ITTK &knownbeg, const ITTK &knownend,
             const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
             const std::function < SCAL (const typename MESHINTERFACE::vertex &) > &Ffunc,
             SCAL epsilon_ratio,
             DATAMANAGER < xgraph::nodeTo < typename MESHINTERFACE::vertex,SCAL,VECT > > &flux_to
             )
{
    DATAMANAGER < xgraph::nodeFrom < typename MESHINTERFACE::vertex,SCAL,VECT,3 > > flux_from;
    internal::FMupdaterFlux < DATAMANAGER, MESHINTERFACE, VECT, SCAL > updater(mi, ls, Ffunc, epsilon_ratio, flux_from, flux_to);
    internal::fmeikFlux(mi, knownbeg, knownend, knownvolatilbeg,knownvolatilend, updater);
}
#endif

} // end of namespace
#endif
