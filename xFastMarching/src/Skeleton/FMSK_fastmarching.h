/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef _FMSK_fastmarching_h_
#define _FMSK_fastmarching_h_
#include <iostream>
#include <array>
#include <cmath>
#include <cassert>
#include <sstream>
#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <limits>
#include <unordered_map>
#include "FM_seq.h"
#include "linearalgebra3d.h"

namespace xfastmarching
{

  namespace skeleton
  {

    template <class SCAL, class VECT2D, class VECT>
      bool  computeTrial2(const VECT2D &T,  const tensor3d2d<SCAL > &dG,  const   SCAL & F, SCAL &Trial, VECT & grad ){
      const VECT2D ones ={1.,1.};
      const VECT dGones = dG*ones;
      const VECT dGT =   dG*T;
      const SCAL A = dot(dGones, dGones);
      const SCAL B = dot(dGones, dGT);
      const SCAL C = dot(dGT, dGT) - F*F;
      SCAL r1;
      if (solvepoly2b(A, B, C, r1) && (r1>=T(0) || r1>=T(1))){
        // if (solvepoly2b(A, B, C, r1) && (r1>=T(0) && r1>=T(1))){
        grad = dG*VECT2D{T(0)-r1, T(1)-r1};
        const VECT2D proj = grad *dG;
        // std::cout << "r1 " << r1 << std::endl;
        if (proj(0) <= 0. && proj(1) <=0. ){
          // std::cout << " OK " << std::endl;
          Trial = r1;
          return true;
        }
      }
      return false;
    }

    template < class ITTK, class ITVK, class ITK, class FMUPDATER>
      void fmeik_intern2 ( const ITTK &knownbeg, const ITTK &knownend
                           ,const ITVK &knownvolatilbeg, const ITVK &knownvolatilend
                           , ITK knownit,
                           FMUPDATER &updater ){
      using meshinterface_t = typename FMUPDATER::meshinterface_t;
      using vertex =  typename meshinterface_t::vertex;
      using edge =    typename meshinterface_t::edge;
      using scal =    typename FMUPDATER::scal;

      const meshinterface_t& mi = updater.mi;
      const scal L = updater.L;
      auto & ls = updater.ls;
      Status<meshinterface_t, vertex > &status=updater.status;

      std::set< std::pair<double, const  vertex *> , compvalpvertex< vertex > > Trial;
      std::list<const edge *> edges;
      //  edges.reserve(50);

      // Put known values in Known status, in  Trial and in knowit
      std::for_each(knownbeg, knownend,
                    [&status, &L, &ls, &Trial, &knownit](const vertex * v)
                    {
                      assert(status(v) == FMStatus::F );
                      knownit=v;
                      ++knownit;
                      status(v) = FMStatus::K;
                      scal Tv = L;
                      ls.get(*v, Tv);
                      Trial.insert( std::make_pair(Tv,  v ) );
                    }
                    );

      // Put  volatil in Trial section with volatil G status only if not already set as known
      std::for_each(knownvolatilbeg, knownvolatilend,[&status, &Trial, &L, &ls](const vertex * v){
          if(status(v) == FMStatus::F)
            {
              status(v) = FMStatus::G;
              scal Tv = L;
              ls.get(*v, Tv);
              Trial.insert( std::make_pair(Tv,  v ) );
            }
        }
        );
      std::list< const vertex * > updated;
      // loop on trial till all known
      while (!Trial.empty()) {
        // std::cout << "Trial Size " << Trial.size() << std::endl;
        auto it = Trial.begin();
        auto vk = it->second;
        Trial.erase(it);
        if(status(vk)==FMStatus::G) {
          knownit=vk;
          ++knownit;
        }
        status(vk) = FMStatus::K;
        edges.clear();
        getEdges(mi, *vk, std::front_inserter (edges));
        for(auto e :edges){
          const vertex *v;
          getothernodes(mi, *vk, *e, v);
          switch (status(v)){
          case FMStatus::F :{
            status(v) = FMStatus::T;
            scal Tv = L;
            ls.set(*v, Tv);
            updater.updatef( *v, *e );
            updated.push_back( v);
            break;
          }
          case FMStatus::T:{
            scal Told = L;
            ls.get(*v, Told);
            if (updater.updatef( *v, *e )){
              Trial.erase(Trial.find ( std::make_pair(Told,   v)));
              updated.push_back(v);
            }
            break;
          }
          case FMStatus::G:{
            scal Told = L;
            ls.get(*v, Told);
            if (updater.updatefg( *v, *e )){
              status(v) = FMStatus::T;
              Trial.erase(Trial.find ( std::make_pair(Told,   v)));
              updated.push_back(v);
            }
            break;
          }
          default:{
            break;
          }
          }
        }
        while(!updated.empty()){
          // std::cout << "updatede " << updated.size() << std::endl;
          const vertex *vk = updated.front();
          updated.pop_front();
          scal Tk = L;
          ls.get(*vk, Tk);
          Trial.insert( std::make_pair(Tk,  vk ) );
          edges.clear();
          getEdges(mi, *vk, std::front_inserter (edges));
          for(auto e :edges) {
            const vertex *v22;
            getothernodes(mi, *vk, *e, v22);
            switch (status(v22)){
            case FMStatus::T:{
              scal Told = L;
              ls.get(*v22, Told);
              if (updater.updatef( *v22, *e )){
                auto tt = Trial.find ( std::make_pair(Told,   v22));
                if (tt != Trial.end())
                  Trial.erase(tt );
                updated.push_back(v22);
              }
              //break;
            }
            default : {break;}
            }
          }
        }
        //loop on sorted trial
      }
    }
    //384
  }//end namespace skeleton
}//end namespace
#endif
