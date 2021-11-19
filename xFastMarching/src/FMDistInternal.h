/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMDISTINTERNAL_H
#define _FMDISTINTERNAL_H

#include "FMDataExchanger.h"

namespace xfastmarching
{
namespace internal
{


template < class MESHINTERFACE, class ITTK, class ITVK, class ITK, class FMUPDATER >
void fmeikDist ( const MESHINTERFACE &mi
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


    size_t nb_to_set = getNbVertex(mi);

    std::set < pair_t, compvalpvertex < vertex,scal > > Trial;
    std::unordered_map < const vertex *, scal > bTrial;
    std::vector < pair_t > new_known;
    new_known.reserve(nb_to_set);
    std::vector < const edge * > edges;
    edges.reserve(50);

    FMDataExchanger < MESHINTERFACE,FMUPDATER > exchange_manager(mi,updater,bTrial);
    MPI_Comm world = exchange_manager.getComm();

#ifdef TIME_MONITORING_FM
    xtool::xDeltaTime dt(world);
    int iddt0 = dt.start("Set K");
    int iddt1 = dt.initAccu("Trial loop");
    int iddt2 = dt.initAccu("Exchange");
    int iddt3 = dt.initAccu("Reset known");
    int iddt4 = dt.initAccu("Reduce");
#endif

    // Put known values in Known status, in  Trial and in knowit
    std::for_each(knownbeg, knownend,
                  [&status, &L, &Trial, &knownit,&updater](const vertex * v)
                  {
                      assert(status(v) == FMStatus::F );
                      knownit = v;
                      ++knownit;
                      status(v) = FMStatus::KF;
                      scal Tv = L;
                      updater.setLsKnown (*v, Tv);
                      Trial.insert( std::make_pair(Tv,  v ) );
                  }
                  );


    /*
       // Put  volatil in Trial section with volatil G status only if not already set as known
       std::for_each(knownvolatilbeg, knownvolatilend,[&status, &Trial, &L, &ls](const vertex * v){
       if (status(v) == FMStatus::F)
       {
       status(v) = FMStatus::G;
       scal Tv = L;
       ls.get(*v, Tv);
       Trial.insert( std::make_pair(Tv,  v ) );
       }
       }
       );
     */

#ifdef TIME_MONITORING_FM
    dt.end(iddt0);
#endif

    // setting for comunication loop
    int loop = 0;
    char not_finish;
    size_t trial_size = Trial.size();
    // 10 here is completelly arbitrarely but avoid first to divide by 0 and hopfully give a not to bad estimate of pass number needed localy
    // to complete FM locally
    //int nb_local_loop = static_cast < int >( nb_to_set/std::max(size_t(10),trial_size));
    int nb_local_loop = ( trial_size > 0 ) ? static_cast < int >( nb_to_set/std::max(size_t(10),trial_size)) : 10;
    MPI_Allreduce(MPI_IN_PLACE,&nb_local_loop,1,MPI_INT,MPI_MAX,world);
    cout<<"nb_local_loop "<<nb_local_loop<<endl;

    //nb_local_loop = 1;
    //cout<<"nb_local_loop "<<nb_local_loop<<endl;

    // loop till Trial empty in all proc
    do
    {

        cout<<"new local loop "<<endl;
        // loop arbitrarely nb_local_loop without reduce
        for (int l = 0; l < nb_local_loop; ++l)
        {
            // reset band to investigate by reseting to zero limit and take actual
            // Trial size: more or less a band depending on update front topology and comunication
            trial_size = Trial.size();
            size_t limit = 0;
            exchange_manager.clearKeyContainer();

#ifdef TIME_MONITORING_FM
            dt.startAccu(iddt1);
#endif

            // loop on trial till empty or limit
            while (!Trial.empty())
            {
                if (++limit > trial_size) break;                  // Band limite to force comunication
                auto it = Trial.begin();
                auto vk = it->second;
                if (status(vk) != FMStatus::KF)
                {
                    new_known.push_back(*it);
                    status(vk) = FMStatus::K;
                }
                Trial.erase(it);
                if (status(vk) == FMStatus::G)
                {
                    knownit = vk;
                    ++knownit;
                }
                // add to exchange
                exchange_manager.accumulate(*vk);
                edges.clear();
                getEdges(mi, *vk, std::back_inserter (edges));
                FM_OUTPUT(" New know "<<vk->getId()<<" with "<<it->first<<endl);
                for (auto e : edges)
                {
                    const vertex *v;
                    getothernodes(mi, *vk, *e, v);
                    switch (status(v))
                    {
                        case FMStatus::F :
                         {
                             FM_OUTPUT(" update Far "<<v->getId()<<endl);
                             status(v) = FMStatus::T;
                             scal Tv;
                             if ( !updater.getGlobalLs(*v, Tv) ) throw -145;
                             updater.setLs(*v, Tv);
                             if ( updater.updatef( *v, *e ) )
                             {
                                 ;
                             }
                             updater.getLs(*v, Tv);
                             Trial.insert( std::make_pair(Tv,  v ) );
                             break;
                         }
                        case FMStatus::T :
                         {
                             FM_OUTPUT(" update Try "<<v->getId()<<endl);
                             scal Told = L;
                             if ( !updater.getGlobalLs(*v, Told) ) throw -145;
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
                             throw -780;
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
                             FM_OUTPUT(" skip "<<v->getId()<<endl);
                             break;
                         }
                    }
                }
                FM_OUTPUT(" End update "<<endl);
            }

#ifdef TIME_MONITORING_FM
            dt.endAccu(iddt1);
            dt.startAccu(iddt2);
#endif

            // Exchange
            scal min_exchanged = exchange_manager.exchange();
            scal o = min_exchanged;
            if (min_exchanged == L) o = 0;
            std::cout.precision(8);
            cout<<"Exchanger "<<++loop<<" "<<scientific<<o<<endl;

#ifdef TIME_MONITORING_FM
            dt.endAccu(iddt2);
#endif
            // Reset known & trial if some significative exchange occures
            if (min_exchanged < L )
            {
#ifdef TIME_MONITORING_FM
                dt.startAccu(iddt3);
#endif
                // container to store new seed if any
                std::set < const vertex * > new_s;
                std::list < const vertex * > bnd;

                FM_OUTPUT(" blunt Trial removal "<<Trial.size()<<endl);
                for (auto &p : Trial)
                {
                    auto vk = p.second;
                    edges.clear();
                    getEdges(mi, *vk, std::back_inserter (edges));
                    for (auto e : edges)
                    {
                        const vertex *v;
                        getothernodes(mi, *vk, *e, v);
                        auto s = status(v);
                        if (s == FMStatus::K || s == FMStatus::KF  )
                            new_s.insert(v);
                    }

                    if (status(vk) != FMStatus::KF  )
                    {
                        updater.setFar(*vk);
                        if (updater.isBoundary(*vk)) bnd.push_back(vk);
                    }
                }
                Trial.clear();


                // clean new_known if any by removing above min_exchanged
                if ( !new_known.empty())
                {
#ifdef OUTPUT_FM
                    for (auto p : new_known)
                    {
                        scal Tv = L;
                        updater.getGlobalLs(*p.second, Tv);
                        cout<<"bbbbb "<<scientific<<p.second->getId()<<" "<<p.first<<" "<<Tv<<endl;

                    }
#endif

                    auto itf = std::find_if( new_known.rbegin(),new_known.rend(),[&min_exchanged,&status,&Trial,&updater,&bnd](pair_t p) -> bool
                                             {
                                                 auto &v = *( p.second );
                                                 if (p.first > min_exchanged)
                                                 {
                                                     assert (status(p.second) == FMStatus::K);
                                                     updater.setFar(v);
                                                     if (updater.isBoundary(v)) bnd.push_back(p.second);

                                                     return false;
                                                 }
                                                 else if (p.first < min_exchanged)
                                                     return true;
                                                 else
                                                 {
                                                     assert (status(p.second) == FMStatus::K);
                                                     updater.setFar(v);
                                                     if (updater.isBoundary(v)) bnd.push_back(p.second);

                                                     return false;
                                                 }
                                             });
                    // all new_know are  greater then min_exchanged so all container must be cleanned
                    if (itf == new_known.rend())
                    {
                        new_known.clear();
                        std::for_each(knownbeg, knownend,
                                      [&status, &new_s](const vertex * v)
                                      {
                                          assert(status(v) == FMStatus::KF );
                                          new_s.insert(v);
                                      }
                                      );
                        assert(!Trial.size());
                    }
                    // all new_know are less or egal to min_exchanged so container must not be changed
                    else if (itf == new_known.rbegin())
                        assert(new_s.size());
                    // remove old known above min_exchanged
                    else
                    {
                        // create range to remove
                        auto itfi = std::next(itf).base();
                        ++itfi;
                        assert(itfi != new_known.end());
                        auto r = xtool::make_range(itfi,new_known.end());

                        // Find Know connected to removed ones
                        for (auto p : r)
                        {
                            auto vk = p.second;
                            edges.clear();
                            getEdges(mi, *vk, std::back_inserter (edges));
                            for (auto e : edges)
                            {
                                const vertex *v;
                                getothernodes(mi, *vk, *e, v);
                                if (status(v) == FMStatus::K || status(v) == FMStatus::KF  )
                                    new_s.insert(v);
                            }
                        }

#ifdef OUTPUT_FM
                        cout<<"to remove "<<r.size()<<" old "<<new_known.size()<<endl;
                        for (auto p : r)
                        {
                            scal Tv = L;
                            updater.getGlobalLs(*p.second, Tv);
                            cout<<"ccccc "<<scientific<<p.second->getId()<<" "<<p.first<<" "<<Tv<<endl;

                        }
#endif
                        new_known.erase(r.begin(),r.end());

                        //
                    }

#ifdef OUTPUT_FM
                    for (auto p : new_known)
                    {
                        scal Tv = L;
                        updater.getGlobalLs(*p.second, Tv);
                        cout<<"aaaaa "<<scientific<<p.second->getId()<<" "<<p.first<<" "<<Tv<<endl;

                    }
#endif
                }

                //Add to trial the boundary received in this pass
                for (auto &p : bTrial)
                {
                    const vertex *v = p.first;
#ifdef OUTPUT_FM
                    {
                        scal Tv = L;
                        updater.getGlobalLs(*v, Tv);
                        cout<<"bTrial "<<scientific<<v->getId()<<" "<<p.second<<" "<<Tv<<" "<<status(v)<<endl;
                    }
#endif
                    switch (status(v))
                    {
                        case FMStatus::T :
                         {
                             scal Tnew,Told = L;
                             if ( !updater.getGlobalLs(*v, Tnew) ) throw -145;
                             auto itf = std::find_if( Trial.begin(),Trial.end(),[&v,&Told](pair_t p) -> bool
                                                      {
                                                          if (p.second == v)
                                                          {
                                                              Told = p.first;
                                                              return true;
                                                          }
                                                          else
                                                              return false;
                                                      });
                             if (itf != Trial.end() && Tnew < Told)
                             {
                                 Trial.erase(itf);
                                 Trial.insert( std::make_pair(Tnew,  v ) );
                                 updater.setLs(*v, Tnew);
                             }
                             break;
                         }
                        case FMStatus::K :
                         {
                             updater.setFar(*v);
                         }
                        case FMStatus::F :
                         {
                             scal Tnew;
                             status(v) = FMStatus::T;
                             if ( !updater.getGlobalLs(*v, Tnew) ) throw -147;
                             updater.setLs(*v, Tnew);
                             Trial.insert( std::make_pair(Tnew,  v ) );
                             edges.clear();
                             getEdges(mi, *v, std::back_inserter (edges));
                             for (auto e : edges)
                             {
                                 const vertex *vo;
                                 getothernodes(mi, *v, *e, vo);
                                 if (status(vo) == FMStatus::K || status(vo) == FMStatus::KF  )
                                     new_s.insert(vo);
                             }
                             break;
                         }
                        default :
                            FM_OUTPUT("bTrial status "<<scientific<<v->getId()<<" "<<status(v)<<endl);
                            throw -148;
                    }
                }
                bTrial.clear();

                // Use seed found above to populate trial
                for (auto vk : new_s)
                {
                    // if seed is still a K we update around
                    if (status(vk) == FMStatus::K || status(vk) == FMStatus::KF  )
                    {
                        edges.clear();
                        getEdges(mi, *vk, std::back_inserter (edges));
                        FM_OUTPUT(" seed know "<<vk->getId()<<endl);
                        for (auto e : edges)
                        {
                            const vertex *v;
                            getothernodes(mi, *vk, *e, v);
                            switch (status(v))
                            {
                                case FMStatus::F :
                                 {
                                     FM_OUTPUT(" update Far "<<v->getId()<<endl);
                                     status(v) = FMStatus::T;
                                     scal Tv;
                                     if ( !updater.getGlobalLs(*v, Tv) ) throw -145;
                                     updater.setLs(*v, Tv);
                                     updater.updatef( *v, *e );
                                     updater.getLs(*v, Tv);
                                     Trial.insert( std::make_pair(Tv,  v ) );
                                     break;
                                 }
                                case FMStatus::T :
                                 {
                                     FM_OUTPUT(" update Try "<<v->getId()<<endl);
                                     scal Told = L;
                                     if ( !updater.getGlobalLs(*v, Told) ) throw -145;
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
                                     throw -780;
                                     break;
                                 }
                                default : {
                                     FM_OUTPUT(" skip "<<v->getId()<<endl);
                                     break;
                                 }
                            }
                        }
                        FM_OUTPUT(" End update "<<endl);

                    }
                    // else
                    // seed is Far: we do not change it
                    // seed becomes a T we do not change it
                }

                // add to Trial if not already the case bnd that have been cleaned and not forcelly
                // in bTrial due to there lack of updating form another proc.
                for (auto vk : bnd)
                {
                    // if F pass it to T
                    if (status(vk) == FMStatus::F )
                    {
                        FM_OUTPUT(" update bnd Far "<<vk->getId()<<endl);
                        assert(updater.isBoundary(*vk));
                        status(vk) = FMStatus::T;
                        scal Tv;
                        if ( !updater.getGlobalLs(*vk, Tv) ) throw -145;
                        updater.setLs(*vk, Tv);
                        Trial.insert( std::make_pair(Tv,  vk ) );
                    }
                }

#ifdef TIME_MONITORING_FM
                dt.endAccu(iddt3);
#endif
            }
        }

#ifdef TIME_MONITORING_FM
        dt.startAccu(iddt4);
#endif
        not_finish = ( Trial.size() > 0 ) ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE,&not_finish,1,MPI_CHAR,MPI_MAX,world);

#ifdef TIME_MONITORING_FM
        dt.endAccu(iddt4);
#endif


    } while (not_finish);

#ifdef TIME_MONITORING_FM
    dt.print();
#endif

}
} // enf of namespace xfastmarching
} // enf of namespace xfastmarching
#endif
