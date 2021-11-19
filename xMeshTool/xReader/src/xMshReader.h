/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */

#ifndef _MSHREADER_XMESHTOOL
#define _MSHREADER_XMESHTOOL
#include<iostream>

namespace xmeshtool
{
  /// A basic msh reader that gives the mesh in theMesh and a association of readed element and their id in the file in entities_id
  /*!
   * This the generic template version
   * Specialisation must be settle in appropriate location (i.e. for now in xinterface::aomd interface for AOMD mesh)
   * M is the type of mesh data base
   * A must follow DATAMANAGER requirements (see xUnorderedMapDataManager.h, xGeneralUnorderedMapDataManager.h and xAttachedDataManagerAOMD.h)
   */
  template<typename M, typename A>
  void xMshReader(std::istream &f,M &theMesh, A & entities_id);

  /// A basic msh reader that gives the mesh in theMesh 
  /*!
   * This the generic template version
   * Specialisation must be settle in appropriate location (i.e. for now in xinterface::aomd interface for AOMD mesh)
   * M is the type of mesh data base
   */
  template<typename M>
  void xMshReader(std::istream &f,M &theMesh);
} // end namespace

#endif

