/*  
    xfem : C++ Finite Element Library 
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/
#include "xAspectRatioAOMD.h"
// STD
#include <cassert>
#include <limits>

// AOMD
#include "mEntity.h"
#include "mEdge.h"
#include "mFace.h"
#include "mTet.h"
#include "mVertex.h"

#define SCALE23 0.81649658092772

using namespace AOMD;


namespace xtool
{
 
 
double xAspectRatio<AOMD::mTet,enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >::aspectRatio (AOMD::mEntity * entity)
{
    assert(dynamic_cast < AOMD::mTet * >( entity ));
    return xAspectRatio<mTet,enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >::aspectRatio(static_cast < AOMD::mTet * >(entity));
}
double xAspectRatio<mTet,enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >::aspectRatio (AOMD::mTet * entity)
{
    assert(entity->size(0)==4);
    return xAspectRatio<mTet,enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >::aspectRatio(entity->get(0,0),entity->get(0,1),entity->get(0,2),entity->get(0,3));
}
double xAspectRatio<mTet,enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >::aspectRatio (AOMD::mEntity *v1, AOMD::mEntity *v2,AOMD::mEntity *v3, AOMD::mEntity *v4)
{
    assert(dynamic_cast < AOMD::mVertex * >( v1 ));
    assert(dynamic_cast < AOMD::mVertex * >( v2 ));
    assert(dynamic_cast < AOMD::mVertex * >( v3 ));
    assert(dynamic_cast < AOMD::mVertex * >( v4 ));
    std::vector < Trellis_Util::mPoint > p;
    p.reserve(4);
    p.push_back( static_cast < AOMD::mVertex * >( v1 )->point());
    p.push_back( static_cast < AOMD::mVertex * >( v2 )->point());
    p.push_back( static_cast < AOMD::mVertex * >( v3 )->point());
    p.push_back( static_cast < AOMD::mVertex * >( v4 )->point());

    int i,i1,i2,i3;
    double v01[3];
    double v02[3];
    double v03[3];
    double n012[3];
    double norm012,norm01,norm02,height;
    double mx=0.;
    double mi=std::numeric_limits<double>::max();

    for (i = 0; i < 4; ++i)
    {
        i1 = ( i+1 )%4;
        i2 = ( i+2 )%4;
        i3 = ( i+3 )%4;
        v01[0] = p[i1](0)-p[i](0);
        v01[1] = p[i1](1)-p[i](1);
        v01[2] = p[i1](2)-p[i](2);
        v02[0] = p[i2](0)-p[i](0);
        v02[1] = p[i2](1)-p[i](1);
        v02[2] = p[i2](2)-p[i](2);
        v03[0] = p[i3](0)-p[i](0);
        v03[1] = p[i3](1)-p[i](1);
        v03[2] = p[i3](2)-p[i](2);
        n012[0] = v01[1]*v02[2]-v01[2]*v02[1];
        n012[1] = v01[2]*v02[0]-v01[0]*v02[2];
        n012[2] = v01[0]*v02[1]-v01[1]*v02[0];
        norm012 = sqrt(n012[0]*n012[0]+n012[1]*n012[1]+n012[2]*n012[2]);
        norm01  = v01[0]*v01[0]+v01[1]*v01[1]+v01[2]*v01[2];
        norm02  = v02[0]*v02[0]+v02[1]*v02[1]+v02[2]*v02[2];

        height = fabs(v03[0]*n012[0]+v03[1]*n012[1]+v03[2]*n012[2])/norm012;
        mi=(height<mi)?height:mi;
        norm01=(norm01>norm02)?norm01:norm02;
        mx=(norm01>mx)?norm01:mx;
    }

    return SCALE23*sqrt(mx)/mi;
}

} // end of namespace



