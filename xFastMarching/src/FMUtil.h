/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMUTIL_H
#define _FMUTIL_H


//#define TIME_MONITORING_FM

#ifdef TIME_MONITORING_FM
#include "xDeltaTime.h"
#endif

//#define OUTPUT_FM
#ifdef OUTPUT_FM
#define FM_OUTPUT(a) cout<<a
#else
#define FM_OUTPUT(a)
#endif

#include "linearalgebra3d.h"

namespace xfastmarching
{

// ==================================================================================
// Class representing Transport concept requirements
template < class SCAL >
class transportConcept
{
    public:
        transportConcept operator*(const SCAL &rhs) const {return *this; }
        transportConcept operator/(const SCAL &rhs) const {return *this; }
        transportConcept & operator/=(const SCAL &scal) {return *this; }
        transportConcept & operator=(const transportConcept &in){return *this; }
        transportConcept operator+(const transportConcept & rhs) const {return *this; }

};

// ==================================================================================

namespace internal
{

enum FMStatus {F,T,K,KF,G};

// solve AX2 -2BX +C
template < class SCAL >
inline bool solvepoly2b(SCAL a, SCAL b, SCAL c, SCAL &r1)
{
    SCAL delta = b*b-a*c;
    if (delta < 0) return false;
    else
    {
        SCAL sqrtd = sqrt(delta);
        r1 = ( b +sqrtd )/a;
        return true;
    }
}



template < class SCAL, class VECT >
bool  computeTrial(const SCAL  &T,  const VECT  &dbase,  const SCAL & F, SCAL &trial, VECT & grad )
{
    trial = T + F/nrm2(dbase);
    grad = dbase*( T - trial );
    return true;
}


template < class SCAL, class VECT2D, class VECT >
bool  computeTrial(const VECT2D &T,  const tensor3d2d < SCAL > &dG,  const SCAL & F, SCAL &Trial, VECT & grad )
{
    const VECT2D ones = {xtool::xDataType < SCAL >::one(),xtool::xDataType < SCAL >::one()};
    const VECT dGones = dG*ones;
    const VECT dGT = dG*T;
    const SCAL A = dot(dGones, dGones);
    const SCAL B = dot(dGones, dGT);
    const SCAL C = dot(dGT, dGT) - F*F;
    const SCAL zero = xtool::xDataType < SCAL >::zero();
    SCAL r1;
    if (solvepoly2b(A, B, C, r1) && r1 >= T(0) && r1 >= T(1))
    {
        grad = dG*VECT2D {T(0)-r1, T(1)-r1};
        const VECT2D proj = grad *dG;
        if (proj(0) <= zero && proj(1) <= zero )
        {
            Trial = r1;
            return true;
        }
    }
    return false;
}

template < class SCAL, class VECT >
bool  computeTrial(const VECT &T,  const tensor3d < SCAL > &dG,  const SCAL & F, SCAL &Trial, VECT & grad)
{
    const VECT ones = {xtool::xDataType < SCAL >::one(),xtool::xDataType < SCAL >::one(),xtool::xDataType < SCAL >::one()};
    const VECT dGones = dG*ones;
    const VECT dGT = dG*T;
    const SCAL A = dot(dGones, dGones);
    const SCAL B = dot(dGones, dGT);
    const SCAL C = dot(dGT, dGT) - F*F;
    const SCAL zero = xtool::xDataType < SCAL >::zero();
    SCAL r1;
    if (solvepoly2b(A, B, C, r1) && ( r1 >= T[0] ) && ( r1 >= T[1] ) && ( r1 >= T[2] ))
    {
        grad = dG*VECT {T[0]-r1, T[1]-r1, T[2]-r1};
        const VECT proj = ( grad * dG );
        if (( proj(0) <= zero  && proj(1) <= zero && proj(2) <= zero ))
        {
            Trial = r1;
            return true;
        }
    }
    return false;
}

template < class MESHINTERFACE >
inline void getothernodes(const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::edge &e,
                          const typename MESHINTERFACE::vertex  * & va)
{
    const typename MESHINTERFACE::vertex *  vv[2];
    getVertices(mi, e, vv );
    if (vv[1] == &v) {va = vv[0]; return; }
    if (vv[0] == &v) {va = vv[1]; return; }
    throw;
}

template < class MESHINTERFACE >
inline const typename MESHINTERFACE::vertex * getothernodes(const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::edge &e)
{
    const typename MESHINTERFACE::vertex *ret = nullptr;
    getothernodes(mi, v, e, ret);
    return ret;
}

template < class MESHINTERFACE >
inline void getothernodes(const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::face &f,
                          const typename MESHINTERFACE::vertex  * & va, const typename MESHINTERFACE::vertex * & vb)
{
    const typename MESHINTERFACE::vertex *  vv[3];
    getVertices(mi, f, vv );
    if (vv[0] == &v) {va = vv[1]; vb = vv[2]; return; }
    if (vv[1] == &v) {va = vv[2]; vb = vv[0]; return; }
    if (vv[2] == &v) {va = vv[0]; vb = vv[1]; return; }
    throw;
}

template < class MESHINTERFACE >
std::array < const typename MESHINTERFACE::vertex *, 2 > getothernodes( const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::face &f)
{
    std::array < const typename MESHINTERFACE::vertex *, 2 > res {nullptr, nullptr };
    getothernodes(mi, v, f, res[0], res[1]);
    return res;
}


template < class MESHINTERFACE >
inline void getothernodes( const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::region &r,
                           const typename MESHINTERFACE::vertex  * & va, const typename MESHINTERFACE::vertex * & vb, const typename MESHINTERFACE::vertex * &vc)
{
    const typename MESHINTERFACE::vertex *  vv[4];
    getVertices(mi, r,  vv );
    if (vv[0] == &v) {va = vv[1]; vb = vv[2]; vc = vv[3];  return; }
    if (vv[1] == &v) {va = vv[2]; vb = vv[3]; vc = vv[0]; return; }
    if (vv[2] == &v) {va = vv[3]; vb = vv[0]; vc = vv[1]; return; }
    if (vv[3] == &v) {va = vv[0]; vb = vv[1]; vc = vv[2]; return; }
    throw;
}

template < class MESHINTERFACE >
std::array < const typename MESHINTERFACE::vertex *, 3 > getothernodes( const MESHINTERFACE &mi, const typename MESHINTERFACE::vertex &v, const typename MESHINTERFACE::region &r)
{
    std::array < const typename MESHINTERFACE::vertex *, 3 > res {nullptr, nullptr, nullptr};
    getothernodes(mi, v, r, res[0], res[1], res[2]);
    return res;
}

template < typename VERTEX, typename SCAL >
struct compvalpvertex
{
    bool operator()( const std::pair < SCAL,  const VERTEX * > & v1, const std::pair < SCAL,  const VERTEX * > & v2 ) const
    {
        if (v1.first < v2.first) return true;
        if (v1.first > v2.first) return false;
        if (v1.second < v2.second) return true;
        return false;
    }
};

template < class MESHINTERFACE, class VERTEX >
class Status
{
    public:
        Status(const MESHINTERFACE &_mi) : status_tag(0), mi(_mi)
        {
            std::stringstream tagname;
            tagname << "STATUS" << this;
            status_tag = getNewTag(mi, tagname.str());
        }
        FMStatus &operator()(const VERTEX *v)
        {
            void *s = getAttachedDataPointer(mi, const_cast < VERTEX & >( *v ), status_tag);
            if (s)
            {
                return *( static_cast < FMStatus * >( s ));
            }
            else
            {
                slist.push_back(FMStatus::F);
                vlist.push_back(v);
                attachDataPointer(mi, const_cast < VERTEX & >( *v ), status_tag, &slist.back() );
                return slist.back();
            }
        }
        FMStatus operator()(const VERTEX *v) const
        {
            void *s = getAttachedDataPointer(mi, const_cast < VERTEX & >( *v ), status_tag);
            if (s) return *( static_cast < FMStatus * >( s ));
            return FMStatus::F;
        }
        void clear()
        {
            std::size_t status_tag_ = status_tag;
            for_each(vlist.begin(), vlist.end(),[&status_tag_](const VERTEX *v) {const_cast < VERTEX * >( v )->deleteData(status_tag_); });
            vlist.clear();
            slist.clear();
        }
        ~Status()
        {
            clear();
            releaseTag(mi, status_tag);
        }

    private:
        std::list < const VERTEX * > vlist;
        std::list < FMStatus > slist;
        std::size_t status_tag;
        const MESHINTERFACE &mi;
};


} // end namespace internal
} // end namespace xfastmarching

#endif
