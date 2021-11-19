/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMDIST_H_
#define _FMDIST_H_

#include "FMUtil.h"
#include "FMDistInternal.h"
#include "FMUpdaterDist.h"

namespace xfastmarching
{


/// Driver for fast marching problem.
/*!
   Compute in ls result of fast marching starting with values at known vertices (given by knownbeg and knownend) ls values at the nodes in this iterator ranges must be known inside ls. "Velocity " or target norm of the gradient at nodes must be given via Ffunc.
   All queries to meshdatabase are done throught the MESHINTERFACE parameter.
   Parameters gls and vn correspond respectively to ls gradient and transported quantity. They are infact kind of optional :  if not given one of the version below is actually launched.
   vn can be represented by any type T defined at nodes, via entity storage. They must be given at "known" nodes. type T must overload left multiplication by double : T (double * T), addition :  T ( T+T), assignement T = T.
   gls is represented as a "vector2d"
   !*/
template < template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE, class ITTK, class ITVK, class ITK, class SCAL, class TRANSPORTED >
void fmeikDist ( const MESHINTERFACE &mi,
                 DATAMANAGER < SCAL > &ls,
                 const ITTK &knownbeg, const ITTK &knownend,
                 const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
                 ITK knownit,
                 const std::function < SCAL (const typename MESHINTERFACE::vertex &) > &Ffunc,
                 SCAL epsilon_ratio,
                 DATAMANAGER < VECT > &gls,
                 DATAMANAGER < TRANSPORTED > &vn)
{
    internal::FMupdaterGenericDist < DATAMANAGER, MESHINTERFACE, VECT, SCAL, TRANSPORTED > updater(mi, ls, Ffunc, epsilon_ratio, gls, vn);
    internal::fmeikDist(mi, knownbeg, knownend, knownvolatilbeg,knownvolatilend,knownit, updater);
}

/// This version does the fast marching but only compute the gls (levelset gradient at node) additional quantity.
template < template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE, class ITTK, class ITVK, class ITK, class SCAL >
void fmeikDist ( const MESHINTERFACE &mi,
             DATAMANAGER < SCAL > &ls,
             const ITTK &knownbeg, const ITTK &knownend,
             const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
             ITK knownit,
             const std::function < SCAL (const typename MESHINTERFACE::vertex &) > &Ffunc,
             SCAL epsilon_ratio,
             DATAMANAGER < VECT > &gls)
{
    internal::FMupdaterGenericDist < DATAMANAGER, MESHINTERFACE, VECT, SCAL, transportConcept < SCAL > > updater(mi, ls, Ffunc, epsilon_ratio, gls);
    internal::fmeikDist(mi, knownbeg, knownend,knownvolatilbeg,knownvolatilend, knownit, updater);
}

/// This version does the fast marching but without computing any additional quantities.
template < template < class >  class DATAMANAGER, class VECT, class MESHINTERFACE,  class ITTK, class ITVK,  class ITK, class SCAL >
void fmeikDist ( const MESHINTERFACE &mi,
             DATAMANAGER < SCAL > &ls,
             const ITTK &knownbeg, const ITTK &knownend,
             const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
             ITK knownit,
             const std::function < SCAL (const typename MESHINTERFACE::vertex &) > &Ffunc,
             SCAL epsilon_ratio)
{
    internal::FMupdaterGenericDist < DATAMANAGER, MESHINTERFACE, VECT, SCAL, transportConcept < SCAL > > updater(mi, ls, Ffunc, epsilon_ratio );
    internal::fmeikDist(mi, knownbeg, knownend, knownvolatilbeg,knownvolatilend, knownit, updater);
}



} // enf of namespace xfastmarching
#endif
