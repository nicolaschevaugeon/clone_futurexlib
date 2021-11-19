/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _X_MESH_TOOL_AOMD_H_
#define _X_MESH_TOOL_AOMD_H_

#include <vector>

#include "mEdge.h"
#include "mEntity.h"
#include "mFace.h"
#include "mHex.h"
#include "mRegion.h"
#include "mTet.h"
#include "xAttachedDataManagerAOMD.h"

// interface aomd
namespace xfem
{
class xMesh;
}

namespace xfem
{
class xMesh;
}

namespace xinterface
{
namespace aomd
{
void modifyAllState(AOMD::mMesh& mmesh);

void modifyAllStateFalse(AOMD::mMesh& mmesh);

/// fill out with a hard copy of in.
/*!
  link between elements  of in and  out are made using attached Entity with  tags xMesh::has_a_copy_tag, and
  xMesh::is_the_copy_of_tag. by defaults, so that the source mesh is not modifyied at all on exit, the attached entity  to get the
  copied element on the source mesh are deleted. Setting the flag clean_tag_on_source_mesh to false permit to keep then
*/
void xCopyMesh(const xfem::xMesh& in, xfem::xMesh& out);

void xCopyMesh(const xfem::xMesh& in, xfem::xMesh& out, xAttachedDataManagerAOMD<AOMD::mEntity*>& is_copied_to,
               xAttachedDataManagerAOMD<AOMD::mEntity*>& is_the_copy_of);

void Simplexify(AOMD::mMesh&);
void DeleteEntity(AOMD::mMesh&, AOMD::mEntity*);
void QuadToTri(AOMD::mMesh&, AOMD::mFace*);
void HexToTet(AOMD::mMesh&, AOMD::mHex*);
void PrismToTet(AOMD::mMesh&, AOMD::mEntity*);
void PyramidToTet(AOMD::mMesh&, AOMD::mEntity*);

void CreateTetWithPrismVertices(AOMD::mMesh&, int N0, int N1, int N2, int N3, int N4, int N5, AOMD::pGEntity);
void CreateTetWithPyramidVertices(AOMD::mMesh&, int N0, int N1, int N2, int N3, int N4, AOMD::pGEntity);
int IndexOfMinAmong(std::vector<int>& iValue, std::vector<int>& index);
void PrintMesh(AOMD::mMesh&);
/*! A crud function which print a mesh using printEntity function. Only topological information are given.
 * If deep argument is true all declared adjacency of all entities are also printed by calling printEntity
 * with its deep parameter to true.
 */
void printMesh(std::ostream& os, xfem::xMesh& mesh, bool deep);

}  // namespace aomd
}  // namespace xinterface

#endif
