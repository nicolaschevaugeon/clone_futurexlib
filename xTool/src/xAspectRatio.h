/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef  _XASPECTRATIO_H
#define  _XASPECTRATIO_H

#include <iostream>


namespace xtool
{

enum enum_xAspectRatio { MAXLENGTHMINHEIGHTRATIO };

//! Generic function object that do nothing : see specialization in other xAspectRatio file
//! T is the type of element considered
//! ART is the type of aspect ratio conputation done for T type element.
//
//! For now only MAXLENGTHMINHEIGHTRATIO criteria is implemented with  AOMD element of tetrahedron kind. (see xAspectRatioAOMD in AOM interfce)
//! MAXLENGTHMINHEIGHTRATIO correspond to the ratio of the maximum edge length
//! divided by the minimum height (considering the 4 heights and 6 length of a tetrahedron)
//
template < typename T, enum_xAspectRatio ART = enum_xAspectRatio::MAXLENGTHMINHEIGHTRATIO >
class xAspectRatio
{
    public:
        static double aspectRatio (T * entity) {std::cout<<"No aspect ratio computation associated to this type"<<std::endl; return 0.; }
};

} // end of namespace

#endif






