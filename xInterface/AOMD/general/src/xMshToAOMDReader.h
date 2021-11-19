/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */

#ifndef _MSHTOAOMDREADER_
#define _MSHTOAOMDREADER_
#include "xMshReader.h"

namespace xinterface
{
namespace aomd
{
template < class T >
class xAttachedDataManagerAOMD;
}
}
namespace AOMD
{
class mMesh;
}

namespace xmeshtool
{
/// A basic msh reader that gives the mesh in AOMD mMesh format and a association of readed element and their id in the file
/*! It a specialzation of xMshReader general template function
 */
template < >
void xMshReader(std::istream &f,AOMD::mMesh &theMesh, xinterface::aomd::xAttachedDataManagerAOMD < int > & entities_id);
/// A basic msh reader that gives the mesh in AOMD mMesh format 
/*! It a specialzation of xMshReader general template function
 */
template < >
void xMshReader(std::istream &f,AOMD::mMesh &theMesh);

} // end namespace

#endif

