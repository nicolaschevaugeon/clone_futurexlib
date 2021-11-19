/*
  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE & LICENSE files for conditions.
*/

#include "InterfaceOctreeToAOMD.h"

#include <algorithm>
#include <functional>

#include "mEdge.h"
#include "mFace.h"
#include "mTet.h"
#include "oExport.h"
#include "oKeyManager.h"
#include "oOctree.h"
#include "xLevelSet.h"
#include "xMesh.h"

// main
using namespace xfem;
using namespace AOMD;
using namespace ::xoctree;

namespace xinterface
{
namespace xoctree
{
void OctreeToAOMDProjection(const oOctree& octree, const oField& lsoct, xMesh* mesh, xLevelSet& ls)
{
   const oKeyManager& key_manager = lsoct.getKeyManager();
   ls.setSupport(mesh, 0.0);
   for (const oKey* key : key_manager)
   {
      mVertex* v = mesh->getMesh().getVertex(key->getId());
      ls(v) = lsoct.getVal(key);
   }
}

void OctreeToAOMDNoHangingProjection(const oOctree& octree, const oField& lsoct, xMesh* mesh, xLevelSet& ls)
{
   const oKeyManager& key_manager = lsoct.getKeyManager();
   // const unsigned int vertex_corner_tag_face =
   //     AOMD_Util::Instance()->lookupMeshDataId("OctreeToAOMDNoHanging_vertex_corner_tag_face");
   //  const unsigned int vertex_corner_tag_vol =
   //     AOMD_Util::Instance()->lookupMeshDataId("OctreeToAOMDNoHanging_vertex_corner_tag_vol");
   xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>> vertex_corner_face_tagger;
   xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>> vertex_corner_vol_tagger;

   ls.setSupport(mesh, 0.0);
   double ls_value;
   for (const oKey* key : key_manager)
   {
      mVertex* v = mesh->getMesh().getVertex(key->getId());
      ls_value = lsoct.getVal(key);
      ls(v) = ls_value;
      // std::vector<AOMD::mEntity*>* vertex_vector = xinterface::aomd::getAttachedEntitiesVector(v, vertex_corner_tag_face);

      if (vertex_corner_face_tagger.getData(*v))
      {
         std::vector<AOMD::mEntity*>* vertex_vector = vertex_corner_face_tagger.getData(*v);
         double mean_contrib = ls_value / 4.;
         std::vector<AOMD::mEntity*>::const_iterator itve = vertex_vector->end();
         for (std::vector<AOMD::mEntity*>::const_iterator itv = vertex_vector->begin(); itv != itve; ++itv)
            ls(static_cast<mVertex*>(*itv)) += mean_contrib;
      }
      //  vertex_vector = xinterface::aomd::getAttachedEntitiesVector(v, vertex_corner_tag_vol);

      if (vertex_corner_vol_tagger.getData(*v))
      {
         std::vector<AOMD::mEntity*>* vertex_vector = vertex_corner_vol_tagger.getData(*v);
         double mean_contrib = ls_value / 8.;
         std::vector<AOMD::mEntity*>::const_iterator itve = vertex_vector->end();
         for (std::vector<AOMD::mEntity*>::const_iterator itv = vertex_vector->begin(); itv != itve; ++itv)
            ls(static_cast<mVertex*>(*itv)) += mean_contrib;
      }
   }
}

void InterfaceAOMDToOctree(const xMesh& mesh, const xLevelSet& ls, oField& lsoct)
{
   const oKeyManager& key_manager = lsoct.getKeyManager();
   for (oKey* key : key_manager)
   {
      mVertex* v = mesh.getMesh().getVertex(key->getId());
      lsoct.setVal(key, ls(v));
   }
}

}  // namespace xoctree

}  // end namespace xinterface
