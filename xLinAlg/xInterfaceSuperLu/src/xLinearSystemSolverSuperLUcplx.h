/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _SOLVERSUPERLUCPLX_H
#define _SOLVERSUPERLUCPLX_H

#ifdef CSUPERLU
namespace superLuCmplx
{
    typedef struct { float r, i; } complex;
} // end namespace
#endif

#ifdef ZSUPERLU
typedef struct { double r, i; } doublecomplex;
#endif

#endif
