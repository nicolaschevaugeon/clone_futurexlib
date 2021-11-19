/*
  octree is a subproject of  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

// -*- C++ -*-

#ifndef _OOCTREE_EXPORT_H__
#define _OOCTREE_EXPORT_H__

#include <string>
#include <vector>

namespace xoctree
{
class oOctree;
class oLevelSet;
class oKeyManager;
class oCellFieldOnFine;

void ExportGMSHAscii(const oLevelSet& l, const oOctree& octree, const std::string& fn, bool simplex = false);
void ExportGMSHAscii(const oCellFieldOnFine& cellfield, const oOctree& octree, const std::string& fn);
void ExportGMSHBinary(const oLevelSet& l, const oOctree& octree, const std::string& fn);

void ExportENSIGHTAscii(const oLevelSet& l, const oOctree& octree, const std::string& fn);
void ExportENSIGHTBinary(const oLevelSet& l, const oOctree& octree, const std::string& fn);

void ExportGMSHAsciiOctreeLevels(const oOctree& octree, const std::string& fn, bool simplex = false);
void ExportGMSHAsciiNodes(const oOctree& octree, const oKeyManager& key_manager, const std::string& fn);
void ExportGMSHMesh(const oOctree& octree, const oKeyManager& k, const std::string& fn);

void ExportGMSHAsciiCurvature(const oLevelSet& l, const oOctree& octree, const std::string& fn, bool simplex = false);
void ExportGMSHAsciiOctreeOnFinestLevel(const oOctree& octree, const std::vector<double>& valVec, const std::string& fn,
                                        bool simplex = false);

}  // namespace xoctree
#endif
