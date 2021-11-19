/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef _FM_h_
#define _FM_h_
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
#include "linearalgebra3d.h"
#include <functional>


namespace xfastmarching
{

  template< class MESHINTERFACE, class ENTITY, class T >
    class entitystorage{
  public:
  entitystorage(const MESHINTERFACE &_mi):datatag(0), mi(_mi){
      std::stringstream tagname;
      tagname << "entitystorage_" << this;
      datatag = getNewTag(mi, tagname.str());
    }

  entitystorage( const entitystorage<MESHINTERFACE, ENTITY,  T > &in):datatag(0), mi(in.mi){
      std::stringstream tagname;
      tagname << "entitystorage_" << this;
      datatag = getNewTag(mi, tagname.str());
      for_each(in.vlist.begin(), in.vlist.end(), [&in, this] (const ENTITY *v){ T val=0.; in.get(*v, val); set(*v, val); }  );
    }

    const std::list<const ENTITY * > & getSupport() const{
      return vlist;
    }

    /// try to return a reference to the data stored at entity v.
    //*! Throw an exception if no data strored at v.
    const T &operator()(const ENTITY &v ) const{
      void *pp  = getAttachedDataPointer(mi,  const_cast< ENTITY &>(v), datatag);
      assert(pp != nullptr);
      return *static_cast<T *>(pp);
    }

    bool get (const ENTITY &v, T &p) const{
      void *pp  = getAttachedDataPointer(mi,  const_cast< ENTITY &>(v), datatag);
      if (pp) {p = *(static_cast<T *> (pp)); return true; };
      return false;
    }

    T * get (const ENTITY &v) const
    {
      void *pp = getAttachedDataPointer(mi,  const_cast < ENTITY & >( v ), datatag);
      return ( static_cast < T * > ( pp ));
    }

    void set (const ENTITY & v, const T &p){
      void * pp  = getAttachedDataPointer(mi, const_cast<ENTITY &>(v), datatag);
      if (!pp){
        pp = new T(p);
        attachDataPointer(mi, const_cast< ENTITY & > (v), datatag, pp);
        vlist.push_back(&v);
      }
      else  *static_cast < T * > (pp) = p;
    }

    ~entitystorage(){
      std::size_t datatag_ = datatag;
      for_each(vlist.begin(), vlist.end(), [&datatag_, this](const ENTITY *v){
          void *pp =  getAttachedDataPointer( mi, const_cast< ENTITY &  >( *v), datatag_);
          assert(pp);
          delete static_cast< T *>(pp);
          deleteData( mi,  const_cast< ENTITY &  >( *v),  datatag_);
        });
      releaseTag(mi, datatag);
    }

  private:
    std::size_t datatag;
    std::list< const ENTITY *> vlist;
    const MESHINTERFACE &mi;
  };



  // solve AX2 -2BX +C
  template<class SCAL >
    inline bool solvepoly2b(SCAL a, SCAL b, SCAL c, SCAL &r1){
    SCAL delta = b*b-a*c;
    if (delta < 0) return false;
    else {
      SCAL sqrtd = sqrt(delta);
      r1 = (b +sqrtd)/a;
      return true;
    }
  }


  template <class SCAL, class VECT>
    bool  computeTrial(const SCAL  &T,  const VECT  &dbase,  const   SCAL & F, SCAL &trial, VECT & grad ){
    trial = T + F/nrm2(dbase);
    grad = dbase*(T - trial);
    return true;
  }


  template <class SCAL, class VECT2D, class VECT>
    bool  computeTrial(const VECT2D &T,  const tensor3d2d<SCAL > &dG,  const   SCAL & F, SCAL &Trial, VECT & grad ){
    const VECT2D ones ={1.,1.};
    const VECT dGones = dG*ones;
    const VECT dGT =   dG*T;
    const SCAL A = dot(dGones, dGones);
    const SCAL B = dot(dGones, dGT);
    const SCAL C = dot(dGT, dGT) - F*F;
    SCAL r1;
    if (solvepoly2b(A, B, C, r1) && r1>=T(0) && r1>=T(1)){
      grad = dG*VECT2D{T(0)-r1, T(1)-r1};
      const VECT2D proj = grad *dG;
      if (proj(0) <= 0. && proj(1) <=0. ){
        Trial = r1;
        return true;
      }
    }
    return false;
  }

  template <class SCAL, class VECT>
    bool  computeTrial(const VECT &T,  const tensor3d<SCAL> &dG,  const   SCAL & F, SCAL &Trial, VECT & grad){
    const VECT ones ={1.,1.,1.};
    const VECT dGones = dG*ones;
    const VECT dGT =   dG*T;
    const SCAL A = dot(dGones, dGones);
    const SCAL B = dot(dGones, dGT);
    const SCAL C = dot(dGT, dGT) - F*F;
    SCAL r1;
    if (solvepoly2b(A, B, C, r1) && (r1>= T[0]) && (r1>=T[1])&& (r1>=T[2])){
      grad = dG*VECT{T[0]-r1, T[1]-r1, T[2]-r1};
      const VECT proj = (grad * dG);
      if(( proj(0) <= 0.  && proj(1) <= 0. && proj(2) <= 0.)){
        Trial = r1;
        return true;
      }
    }
    return false;
  }

  template <class MESHINTERFACE >
    inline void getothernodes(const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::edge &e,
                              const typename MESHINTERFACE::vertex  * & va){
    const typename MESHINTERFACE::vertex *  vv[2];
    getVertices(mi, e, vv );
    if (vv[1] == &v) {va = vv[0]; return;}
    if (vv[0] == &v) {va = vv[1]; return;}
    throw;
  }

  template <class MESHINTERFACE >
    inline const typename MESHINTERFACE::vertex * getothernodes(const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::edge &e){
    const typename MESHINTERFACE::vertex *ret =nullptr;
    getothernodes(mi, v, e, ret);
    return ret;
  }

  template <class MESHINTERFACE >
    inline void getothernodes(const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::face &f,
                              const typename MESHINTERFACE::vertex  * & va, const typename MESHINTERFACE::vertex * & vb){
    const typename MESHINTERFACE::vertex *  vv[3];
    getVertices(mi, f, vv );
    if (vv[0] == &v) {va = vv[1]; vb = vv[2]; return;}
    if (vv[1] == &v) {va = vv[2]; vb = vv[0]; return;}
    if (vv[2] == &v) {va = vv[0]; vb = vv[1]; return;}
    throw;
  }

  template <class MESHINTERFACE >
    std::array< const typename MESHINTERFACE::vertex *, 2> getothernodes( const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::face &f){
    std::array< const typename MESHINTERFACE::vertex *, 2> res{nullptr, nullptr } ;
    getothernodes(mi, v, f, res[0], res[1]);
    return res;
  }


  template <class MESHINTERFACE >
    inline void getothernodes( const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::region &r,
                               const typename MESHINTERFACE::vertex  * & va, const typename MESHINTERFACE::vertex * & vb, const typename MESHINTERFACE::vertex *&vc){
    const typename MESHINTERFACE::vertex *  vv[4];
    getVertices(mi, r,  vv );
    if (vv[0] == &v) {va = vv[1]; vb = vv[2]; vc = vv[3];  return;}
    if (vv[1] == &v) {va = vv[2]; vb = vv[3]; vc = vv[0]; return;}
    if (vv[2] == &v) {va = vv[3]; vb = vv[0]; vc = vv[1]; return;}
    if (vv[3] == &v) {va = vv[0]; vb = vv[1]; vc = vv[2]; return;}
    throw;
  }

  template <class MESHINTERFACE >
    std::array< const typename MESHINTERFACE::vertex *, 3> getothernodes( const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::region &r){
    std::array< const typename MESHINTERFACE::vertex *, 3> res{nullptr, nullptr, nullptr} ;
    getothernodes(mi, v, r, res[0], res[1], res[2]);
    return res;
  }


  template <class VERTEX>
    struct compvalpvertex{
      bool operator()( const std::pair<double,  const VERTEX * > & v1, const std::pair<double,  const VERTEX * > & v2 ) const{
        if (v1.first < v2.first) return true;
        if (v1.first > v2.first) return false;
        if (v1.second < v2.second) return true;
        return false;
      }
    };

  enum FMStatus{F,T,K,G};

  template <class MESHINTERFACE, class VERTEX >
    class Status{
  public:
  Status(const MESHINTERFACE &_mi):status_tag(0), mi(_mi){
      std::stringstream tagname;
      tagname << "STATUS" << this;
      status_tag = getNewTag(mi, tagname.str());
    }
    FMStatus &operator()(const VERTEX *v){
      void *s = getAttachedDataPointer(mi, const_cast<VERTEX & >(*v), status_tag);
      if (s) {
        return *( static_cast <FMStatus * >(s));
      }
      else {
        slist.push_back(FMStatus::F);
        vlist.push_back(v);
        attachDataPointer(mi, const_cast<VERTEX & >(*v), status_tag, &slist.back() );
        return slist.back();
      }
    }
    void clear(){
        std::size_t status_tag_ =  status_tag;
        for_each(vlist.begin(), vlist.end(),[&status_tag_](const VERTEX *v) {const_cast< VERTEX *>(v)->deleteData(status_tag_);});
        vlist.clear();
        slist.clear();
    }
    FMStatus operator()(const VERTEX *v) const{
      void *s = getAttachedDataPointer(mi, const_cast<VERTEX & >(*v), status_tag);
      if (s) return  *(static_cast <FMStatus * >(s));
      return FMStatus::F;
    }
    ~Status(){
      clear();
      releaseTag(mi, status_tag);
    }

  private:
    std::list<const VERTEX *> vlist;
    std::list<FMStatus> slist;
    std::size_t status_tag;
    const MESHINTERFACE &mi;
  };

    // MESHINTERFACE interface represent the type of the interface to mesh
    // T represent the Type of data that can be transported along the gradient
    // GEOMVECT    is the type used to represent vectors
    template <class MESHINTERFACE, class VECT, class SCAL, class TRANSPORTED >
      class FMupdaterGeneric{
    public:
      typedef SCAL scal;
      typedef VECT vect;
      typedef TRANSPORTED transported;
      typedef MESHINTERFACE meshinterface_t;

      typedef typename MESHINTERFACE::vertex vertex;
      typedef typename MESHINTERFACE::edge edge;
      typedef typename MESHINTERFACE::face face;
      typedef typename MESHINTERFACE::region region;

      typedef std::vector<const region * >  r_container_t;
      typedef std::vector<const face * >    f_container_t;
      typedef entitystorage<  MESHINTERFACE, vertex,  SCAL >        entitystorage_t;
      typedef entitystorage<  MESHINTERFACE, vertex,  vect >        gls_t;
      typedef entitystorage<  MESHINTERFACE, vertex,  TRANSPORTED > vn_t;

          FMupdaterGeneric(const MESHINTERFACE &_mi,  entitystorage_t &_ls,
                           std::function< scal (const vertex &) > _Ffunc, scal _epsilon_ratio):
            mi(_mi), status(mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(nullptr), pvn(nullptr){}
          FMupdaterGeneric(const MESHINTERFACE &_mi,  entitystorage_t &_ls,
                           std::function< scal (const vertex &) > _Ffunc, scal _epsilon_ratio, gls_t &_gls  ):
            mi(_mi), status(mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(&_gls), pvn(nullptr){}
          FMupdaterGeneric(const MESHINTERFACE &_mi,  entitystorage_t &_ls,
                           std::function< scal (const vertex &) > _Ffunc, scal _epsilon_ratio, gls_t &_gls, vn_t &_vn  ):
            mi(_mi), status(mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(&_gls), pvn(&_vn){}

      bool updatefg(const vertex &vx, const edge & e) const {
        const scal Fx = Ffunc(vx);
        const vertex *vb =  getothernodes(mi, vx, e);
        const vect b = getcoord( *vb) -getcoord(vx);
        scal epsilon=0.;
        vect Grad;
        computeTrial(scal(0.), dual(b), Fx, epsilon, Grad);
        epsilon*=epsilon_ratio;
        return updatefp(vx,e,epsilon);
      }

      bool updatef(const vertex &vx, const edge & e) const{
        return updatefp(vx,e,0.);
      }

      const scal L = numeric_limits<scal>::max();
      const MESHINTERFACE & mi;
      Status<MESHINTERFACE, vertex > status;
      entitystorage<  MESHINTERFACE, vertex,  scal > &ls;

      protected:

      std::function< scal (const vertex & )> Ffunc;
      const scal epsilon_ratio;

      gls_t   *pgls;
      vn_t    *pvn;

      virtual bool updatefp(const vertex &vx, const edge & e, const scal &epsilon)  const {
        scal valx = ls( vx);
        vect px = getcoord(vx);
        scal Fx = Ffunc(vx);
        scal Trial;
        vect GradTrial;
        vect Grad;
        TRANSPORTED Trans;
        bool updated = false;
        const std::function<bool ( const vertex *v) > is_known = [this]( const vertex *v) {return status(v) ;};
        std::function<bool (double, double)> less=[epsilon](double v, double w){ return v<w+epsilon; };

        r_container_t  e2r;
        f_container_t  e2f;
        e2r.reserve(50);
        e2f.reserve(50);

        getRegions(mi, e, std::back_inserter(e2r));
        for (auto r : e2r){
          const auto va  = getothernodes(mi, vx, *r);
          if( std::all_of(va.begin(), va.end(), is_known ) ) {
            const vect T = { ls( *va[0]), ls( *va[1]), ls( *va[2])};
            if ( less(T(0), valx) || less(T(1), valx) || less(T(2), valx)){
              const std::array<vect, 3> a ={getcoord(*va[0]) -px, getcoord(*va[1]) -px, getcoord(*va[2])-px  };
              const tensor3d<scal> dualG = dual( tensor3d<scal> {a[0], a[1], a[2]} );
              if (computeTrial(T, dualG, Fx, Trial, GradTrial)){
                if(less(Trial,valx)) {
                  updated =true;
                  valx = Trial;
                  Grad = GradTrial;
                  if (pvn){
                    const std::array<TRANSPORTED, 3> trans_other { (*pvn)(*va[0]), (*pvn)(*va[1]), (*pvn)(*va[2])  };
                    auto gc =  Grad * dualG;
                    Trans = (gc(0) *trans_other[0] + gc(1)*trans_other[1]+gc(2)*trans_other[2])/(gc(0) + gc(1) + gc(2));
                  }
                }
              }
            }
          }
        }

        getFaces(mi,   e, std::back_inserter(e2f));
        for (auto f : e2f){
          const auto va = getothernodes(mi, vx, *f);
          if( std::all_of(va.begin(), va.end(), is_known ) ) {
            const vector2d<scal > T = { ls( *va[0]), ls( *va[1]) };
            if ( less(T(0),valx) || less(T(1), valx) ){
              const std::array<vect, 2> a  {getcoord(*va[0])-px, getcoord(*va[1])-px  };
              const auto  dualG  =  dual (tensor3d2d<scal>( a[0], a[1]));
              if (computeTrial(T, dualG, Fx, Trial, GradTrial)){
                if(less(Trial,valx)) {
                  updated =true;
                  valx = Trial;
                  Grad = GradTrial;
                  if (pvn){
                    const std::array<TRANSPORTED, 2> trans_other { (*pvn)(*va[0]) , (*pvn)(*va[1]) };
                    auto gc =  Grad * dualG;
                    Trans = (gc(0) *trans_other[0] + gc(1)*trans_other[1])/(gc(0)+gc(1));
                  }
                }
              }
            }
          }
        }

        const vertex * vb = getothernodes(mi, vx, e);
        if( is_known(vb)){
          const scal T = ls( *vb);
          if ( less(T, valx) ){
            const vect b = getcoord( *vb) -px;
            if (computeTrial(T, dual(b), Fx, Trial, GradTrial)){
              if(less(Trial,valx)) {
                updated =true;
                valx = Trial;
                Grad = GradTrial;
                if (pvn){
                  Trans = (*pvn)(*vb);
                }
              }
            }
          }
        }

        if (updated){
          ls.set(vx, valx);
          if (pgls) pgls->set(vx, Grad);
          if (pvn)  pvn->set(vx, Trans);
          return true;
        }
        return false;
      }
    public:
      vect getcoord( const vertex & v ) const {
        vect res;
        getCoord(mi, v, res);
        return res;
      }




    };

  template <class MESHINTERFACE, class ITTK, class ITVK, class ITK, class FMUPDATER>
    void fmeik_intern ( const MESHINTERFACE& mi, entitystorage<  MESHINTERFACE, typename MESHINTERFACE::vertex, typename FMUPDATER::scal > & ls
                        ,const ITTK &knownbeg, const ITTK &knownend
                        ,const ITVK &knownvolatilbeg, const ITVK &knownvolatilend
                        , ITK knownit,
                        FMUPDATER &updater ){
    typedef typename MESHINTERFACE::vertex vertex;
    typedef typename MESHINTERFACE::edge edge;
    typedef typename FMUPDATER::scal scal;
    const scal L = updater.L;
    Status<MESHINTERFACE, vertex > &status=updater.status;

    std::set< std::pair<double, const  vertex *> , compvalpvertex< vertex > > Trial;
    std::vector<const edge *> edges;
    edges.reserve(50);

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

    // loop on trial till all known
    while (!Trial.empty()){
      auto it = Trial.begin();
      auto vk = it->second;
      Trial.erase(it);
      if(status(vk)==FMStatus::G) {
        knownit=vk;
        ++knownit;
      }
      status(vk) = FMStatus::K;
      edges.clear();
      getEdges(mi, *vk, std::back_inserter (edges));
      for(auto e :edges){
        const vertex *v;
        getothernodes(mi, *vk, *e, v);
        switch (status(v)){
        case FMStatus::F :{
          status(v) = FMStatus::T;
          scal Tv = L;
          ls.set(*v, Tv);
          updater.updatef( *v, *e );
          ls.get(*v, Tv);
          Trial.insert( std::make_pair(Tv,  v ) );
          break;
        }
        case FMStatus::T:{
          scal Told = L;
          ls.get(*v, Told);
          if (updater.updatef( *v, *e )){
            scal Tnew = L;
            ls.get(*v, Tnew);
            Trial.erase(Trial.find ( std::make_pair(Told,   v)));
            Trial.insert( std::make_pair(Tnew,  v ) );
          }
          break;
        }
        case FMStatus::G:{
          scal Told = L;
          ls.get(*v, Told);
          if (updater.updatefg( *v, *e )){
            status(v) = FMStatus::T;
            scal Tnew = L;
            ls.get(*v, Tnew);
            Trial.erase(Trial.find ( std::make_pair(Told,   v)));
            Trial.insert( std::make_pair(Tnew,  v ) );
          }
          break;
        }
        default:{
          break;
        }
        }
      }
    }
  }
  //384

  /// Driver for fast marching problem.
  /*!
    Compute in ls result of fast marching starting with values at known vertices (given by knownbeg and knownend) ls values at the nodes in this iterator ranges must be known inside ls. "Velocity " or target norm of the gradient at nodes must be given via Ffunc.
    All queries to meshdatabase are done throught the MESHINTERFACE parameter.
    Parameters gls and vn correspond respectively to ls gradient and transported quantity. They are infact kind of optional :  if not given one of the version below is actually launched.
    vn can be represented by any type T defined at nodes, via entity storage. They must be given at "known" nodes. type T must overload left multiplication by double : T (double * T), addition :  T ( T+T), assignement T = T.
    gls is represented as a "vector2d"
    !*/
  template <class VECT, class MESHINTERFACE, class ITTK, class ITVK, class ITK, class SCAL, class TRANSPORTED>
    void fmeik ( const MESHINTERFACE& mi,
                 entitystorage<MESHINTERFACE, typename MESHINTERFACE::vertex, SCAL>& ls,
                 const ITTK &knownbeg, const ITTK &knownend,
                 const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
                 ITK knownit,
                 const std::function<SCAL (const typename MESHINTERFACE::vertex &)>& Ffunc,
                 SCAL epsilon_ratio,
                 entitystorage<MESHINTERFACE, typename MESHINTERFACE::vertex, VECT> & gls,
                 entitystorage<MESHINTERFACE, typename MESHINTERFACE::vertex, TRANSPORTED> &vn){
    FMupdaterGeneric<MESHINTERFACE, VECT, SCAL, TRANSPORTED> updater(mi, ls, Ffunc, epsilon_ratio, gls, vn);
    fmeik_intern(mi, ls, knownbeg, knownend, knownvolatilbeg,knownvolatilend,knownit, updater);
  }

  /// This version does the fast marching but only compute the gls (levelset gradient at node) additional quantity.
  template <class VECT, class MESHINTERFACE, class ITTK, class ITVK, class ITK, class SCAL >
    void fmeik ( const MESHINTERFACE& mi,
                 entitystorage<MESHINTERFACE, typename MESHINTERFACE::vertex, SCAL>& ls,
                 const ITTK &knownbeg, const ITTK &knownend,
                 const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
                 ITK knownit,
                 const std::function<SCAL (const typename MESHINTERFACE::vertex &)>& Ffunc,
                 SCAL epsilon_ratio,
                 entitystorage<MESHINTERFACE, typename MESHINTERFACE::vertex, VECT>& gls){
    FMupdaterGeneric<MESHINTERFACE, VECT, SCAL, int> updater(mi, ls, Ffunc, epsilon_ratio, gls);
    fmeik_intern(mi, ls, knownbeg, knownend,knownvolatilbeg,knownvolatilend, knownit, updater);
  }

  // /// This version does the fast marching but without computing any additional quantities.
  template <class VECT, class MESHINTERFACE,  class ITTK, class ITVK,  class ITK, class SCAL >
    void fmeik ( const MESHINTERFACE& mi,
                 entitystorage<MESHINTERFACE, typename MESHINTERFACE::vertex, SCAL>& ls,
                 const ITTK &knownbeg, const ITTK &knownend, const ITVK &knownvolatilbeg, const ITVK &knownvolatilend,
                 ITK knownit,
                 const std::function<SCAL (const typename MESHINTERFACE::vertex &)>& Ffunc,
                 SCAL epsilon_ratio){
    FMupdaterGeneric<MESHINTERFACE, VECT, SCAL, int> updater(mi, ls, Ffunc, epsilon_ratio );
    fmeik_intern(mi, ls, knownbeg, knownend, knownvolatilbeg,knownvolatilend, knownit, updater);
  }


} // end of namespace
#endif
