/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMINTERNAL_H_
#define _FMINTERNAL_H_


namespace xfastmarching
{
namespace internal
{
template < class MESHINTERFACE, class ITTK, class ITVK, class FMUPDATER >
void fmeikFlux ( const MESHINTERFACE & mi
                 ,const ITTK &knownbeg, const ITTK &knownend
                 ,const ITVK &knownvolatilbeg, const ITVK &knownvolatilend
                 ,FMUPDATER &updater )
{
    typedef typename MESHINTERFACE::vertex vertex;
    typedef typename MESHINTERFACE::edge edge;
    typedef typename FMUPDATER::scal scal;
    typedef std::pair < scal, const vertex * > pair_t;
    const scal L = updater.L;
    Status < MESHINTERFACE, vertex > &status = updater.status;

    std::set < pair_t, compvalpvertex < vertex,scal > > Trial;
    std::vector < const edge * > edges;
    edges.reserve(50);

    // Put known values in  Trial and trial status
    std::for_each(knownbeg, knownend,
                  [&status, &L, &updater, &Trial ](const vertex * v)
                  {
                      assert(status(v) == FMStatus::F );
                      scal Tv = L;
                      updater.getLs(*v, Tv);
                      FM_OUTPUT("initial known "<<v->getId()<<" is set to "<<Tv<<endl);
                      Trial.insert( std::make_pair(Tv,  v ) );
                      status(v) = FMStatus::T;
                  }
                  );

    // Put  volatil in Trial and trial status
    std::for_each(knownvolatilbeg, knownvolatilend,[&status, &Trial, &L, &updater](const vertex * v){
                      if (status(v) == FMStatus::F)
                      {
                          scal Tv = L;
                          updater.getLs(*v, Tv);
                          Trial.insert( std::make_pair(Tv,  v ) );
                          status(v) = FMStatus::T;
                      }
                  }
                  );

    // loop on trial till all known
    while (!Trial.empty())
    {
        auto it = Trial.begin();
        auto vk = it->second;
        Trial.erase(it);
        updater.setToKnown(vk);
        edges.clear();
        getEdges(mi, *vk, std::back_inserter (edges));
        for (auto e : edges)
        {
            const vertex *v;
            getothernodes(mi, *vk, *e, v);
            switch (status(v))
            {
                case FMStatus::F : {
                     status(v) = FMStatus::T;
                     scal Tv = L;
                     updater.setLs(*v, Tv);
                     updater.updatef( *v, *e );
                     updater.getLs(*v, Tv);
                     Trial.insert( std::make_pair(Tv,  v ) );
                     break;
                 }
                case FMStatus::T : {
                     scal Told = L;
                     updater.getLs(*v, Told);
                     if (updater.updatef( *v, *e ))
                     {
                         scal Tnew = L;
                         updater.getLs(*v, Tnew);
                         Trial.erase(Trial.find ( std::make_pair(Told,   v)));
                         Trial.insert( std::make_pair(Tnew,  v ) );
                     }
                     break;
                 }
                default : {
                     break;
                 }
            }
        }
    }
}
template < class MESHINTERFACE, class ITTK, class ITVK, class ITK, class FMUPDATER >
void fmeik ( const MESHINTERFACE & mi
             ,const ITTK &knownbeg, const ITTK &knownend
             ,const ITVK &knownvolatilbeg, const ITVK &knownvolatilend
             , ITK knownit,
             FMUPDATER &updater )
{
    typedef typename MESHINTERFACE::vertex vertex;
    typedef typename MESHINTERFACE::edge edge;
    typedef typename FMUPDATER::scal scal;
    typedef std::pair < scal, const vertex * > pair_t;
    const scal L = updater.L;
    Status < MESHINTERFACE, vertex > &status = updater.status;

    std::set < pair_t, compvalpvertex < vertex,scal > > Trial;
    std::vector < const edge * > edges;
    edges.reserve(50);

    // Put known values in Known status, in  Trial and in knowit
    std::for_each(knownbeg, knownend,
                  [&status, &L, &updater, &Trial, &knownit](const vertex * v)
                  {
                      assert(status(v) == FMStatus::F );
                      knownit = v;
                      ++knownit;
                      status(v) = FMStatus::K;
                      scal Tv = L;
                      updater.getLs(*v, Tv);
                      Trial.insert( std::make_pair(Tv,  v ) );
                  }
                  );

    // Put  volatil in Trial section with volatil G status only if not already set as known
    std::for_each(knownvolatilbeg, knownvolatilend,[&status, &Trial, &L, &updater](const vertex * v){
                      if (status(v) == FMStatus::F)
                      {
                          status(v) = FMStatus::G;
                          scal Tv = L;
                          updater.getLs(*v, Tv);
                          Trial.insert( std::make_pair(Tv,  v ) );
                      }
                  }
                  );

    // loop on trial till all known
    while (!Trial.empty())
    {
        auto it = Trial.begin();
        auto vk = it->second;
        Trial.erase(it);
        if (status(vk) == FMStatus::G)
        {
            knownit = vk;
            ++knownit;
        }
        status(vk) = FMStatus::K;
        edges.clear();
        getEdges(mi, *vk, std::back_inserter (edges));
        for (auto e : edges)
        {
            const vertex *v;
            getothernodes(mi, *vk, *e, v);
            switch (status(v))
            {
                case FMStatus::F : {
                     status(v) = FMStatus::T;
                     scal Tv = L;
                     updater.setLs(*v, Tv);
                     updater.updatef( *v, *e );
                     updater.getLs(*v, Tv);
                     Trial.insert( std::make_pair(Tv,  v ) );
                     break;
                 }
                case FMStatus::T : {
                     scal Told = L;
                     updater.getLs(*v, Told);
                     if (updater.updatef( *v, *e ))
                     {
                         scal Tnew = L;
                         updater.getLs(*v, Tnew);
                         Trial.erase(Trial.find ( std::make_pair(Told,   v)));
                         Trial.insert( std::make_pair(Tnew,  v ) );
                     }
                     break;
                 }
                case FMStatus::G : {
                     scal Told = L;
                     updater.getLs(*v, Told);
                     if (updater.updatefg( *v, *e ))
                     {
                         status(v) = FMStatus::T;
                         scal Tnew = L;
                         updater.getLs(*v, Tnew);
                         Trial.erase(Trial.find ( std::make_pair(Told,   v)));
                         Trial.insert( std::make_pair(Tnew,  v ) );
                     }
                     break;
                 }
                default : {
                     break;
                 }
            }
        }
    }
}

} // end of namespace
} // end of namespace
#endif
