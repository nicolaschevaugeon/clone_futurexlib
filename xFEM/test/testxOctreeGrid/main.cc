#include <fstream>

#include "xMPIEnv.h"
#include "xMappingBuilderHolder.h"
#include "xMesh.h"
#include "xOctreeGrid.h"

using namespace xfem;
using namespace xtensor;
using namespace xtool;
using namespace xgeom;

using myOctantBucket = octantBucket<2, 8, 2>;
using elem_data = myOctantBucket::elem_data;

void exportGMSH(const xPoint &p, std::ostream &out) { out << p[0] << ",  " << p[1] << ", " << p[2]; }

void exportAsCubeGMSH(const xBoundingBox &bb, std::ostream &out)
{
   const auto &min = bb.min;
   const auto &max = bb.max;
   auto dir = max - min;
   std::array<xPoint, 8> cubevert;
   for (int k = 0; k < 2; ++k)
   {
      cubevert[k * 4] = xPoint{min[0], min[1], min[2] + dir[2] * k};
      cubevert[k * 4 + 1] = xPoint{min[0] + dir[0], min[1], min[2] + dir[2] * k};
      cubevert[k * 4 + 2] = xPoint{min[0] + dir[0], min[1] + dir[1], min[2] + dir[2] * k};
      cubevert[k * 4 + 3] = xPoint{min[0], min[1] + dir[1], min[2] + dir[2] * k};
   }
   out << "SH(";
   for (int i = 0; i < 7; ++i)
   {
      exportGMSH(cubevert[i], out);
      out << ", ";
   }
   exportGMSH(cubevert[7], out);
   out << ") {0,0,0,0,0,0,0,0};\n";
}

template <class OCTANT>
void exportGMSH(const OCTANT &buck, std::ostream &out)
{
   out << "View \" octree_grid \" { \n";
   buck.visitLeaves([&out](const OCTANT &b) { exportAsCubeGMSH(b.getBoundingBox(), out); });
   out << "};";
}

int main(int argc, char *argv[])
{
   xMPIEnv::init(argc, argv);
   {  // a hand made case, with 6 "elem", with all there centroid at the same position, to force the octree to go to high level !!
      myOctantBucket buck(xBoundingBox{xPoint{-1., -1., -1.}, xPoint{1., 1., 1.}});
      int elem1, elem2, elem3, elem4, elem5, elem6;
      elem_data data_elem1{&elem1, xPoint{0.1, 0.1, 0.}, xBoundingBox{xPoint{-0.5, -0.5, -0.5}, xPoint{0.5, 0.5, 0.5}}};
      elem_data data_elem2{&elem2, xPoint{0.1, 0.1, 0.}, xBoundingBox{xPoint{-0.5, -0.5, -0.5}, xPoint{0.5, 0.5, 0.5}}};
      elem_data data_elem3{&elem3, xPoint{0.1, 0.1, 0.}, xBoundingBox{xPoint{-0.5, -0.5, -0.5}, xPoint{0.5, 0.5, 0.5}}};
      elem_data data_elem4{&elem4, xPoint{0.1, 0.1, 0.}, xBoundingBox{xPoint{-0.5, -0.5, -0.5}, xPoint{0.5, 0.5, 0.5}}};
      elem_data data_elem5{&elem5, xPoint{0.1, 0.1, 0.}, xBoundingBox{xPoint{-0.5, -0.5, -0.5}, xPoint{0.5, 0.5, 0.5}}};
      elem_data data_elem6{&elem6, xPoint{0.1, 0.1, 0.}, xBoundingBox{xPoint{-0.5, -0.5, -0.5}, xPoint{0.5, 0.5, 0.5}}};
      buck.addElement(data_elem1);
      buck.addElement(data_elem2);
      buck.addElement(data_elem3);
      buck.addElement(data_elem4);
      buck.addElement(data_elem5);
      buck.addElement(data_elem6);
      std::ofstream out("octree_simple.pos");
      exportGMSH(buck, out);
   }
   {
      // A case build out of an xmesh, but the octree is build directly in the code.
      std::cout << "Testing xOctreeGrid construction and exploration" << std::endl;
      xMesh mesh("data/CubeWithCubeHole.msh");
      std::cout << "Mesh read : " << mesh.size(0) << " vertices " << mesh.size(3) << " solids " << std::endl;
      myOctantBucket buck(xBoundingBox{xPoint{-0.03, -0.03, -0.03}, xPoint{10.3, 10.3, 10.3}});
      for (auto preg : mesh.range(3))
      {
         std::unique_ptr<xmapping::xMapping> pmapping(xMappingBuilderHolderSingleton::instance().buildMapping(*preg));
         elem_data edata{preg, pmapping->eval(pmapping->COG()), pmapping->boundingBox()};
         buck.addElement(edata);
      }
      std::ofstream out("octree_3d_from_xmesh.pos");
      exportGMSH(buck, out);
      size_t max_nelement = 0;
      buck.visitLeaves([&max_nelement](const myOctantBucket &b) {
         size_t nb = std::distance(b.getDatas().begin(), b.getDatas().end());
         max_nelement = std::max(max_nelement, nb);
      });
      std::cout << "max element in a octree Leaf " << max_nelement << std::endl;
      std::array<xPoint, 2> targets{xPoint{0.05, 0.05, 0.05}, xPoint{0., 0., 0.}};
      for (auto target : targets)
      {
         const myOctantBucket *pbuck = buck.findBucket(target);
         std::list<std::pair<AOMD::mEntity *, xPoint>> pent_to_loccoord;
         if (pbuck)
         {
            std::cout << " found !" << std::distance(pbuck->getDatas().begin(), pbuck->getDatas().end()) << std::endl;
            for (const auto &edata : pbuck->getDatas())
            {
               if (edata.bb.contains(target))
               {
                  AOMD::mEntity *pe = static_cast<AOMD::mEntity *>(edata.elem);
                  std::unique_ptr<xmapping::xMapping> pmapping(xMappingBuilderHolderSingleton::instance().buildMapping(*pe));
                  xPoint uvw;
                  if (pmapping->interiorCheck(target, uvw))
                  {
                     pent_to_loccoord.push_back(std::make_pair(pe, uvw));
                  }
               }
            }
         }
         std::cout << "Point " << target << " Found in " << pent_to_loccoord.size() << " Elements " << std::endl;
         for (const auto &pe_pt : pent_to_loccoord)
         {
            const AOMD::mEntity *pe = pe_pt.first;
            xPoint uvw = pe_pt.second;
            std::unique_ptr<xmapping::xMapping> pmapping(xMappingBuilderHolderSingleton::instance().buildMapping(*pe));
            xPoint xyz = pmapping->eval(uvw);
            std::cout << "Entity " << pe_pt.first << " local Coordinate " << pe_pt.second << " Global Coordinate " << xyz
                      << std::endl;
         }
      }
   }
   {  // A case on a 3d mesh, using the xMesh internal createOctree and locateOctree
      xMesh mesh("data/CubeWithCubeHole.msh");
      std::cout << "Mesh read : " << mesh.size(0) << " vertices " << mesh.size(3) << " solids " << std::endl;
      std::array<xPoint, 2> targets{xPoint{0.05, 0.05, 0.05}, xPoint{0., 0., 0.}};
      std::cout << "Testing mesh.locateElementOctree" << std::endl;
      for (auto target : targets)
      {
         auto entityfounds = mesh.locateElementOctree(target);
         std::cout << "Point " << target << " Found in " << entityfounds.size() << " Elements " << std::endl;
         for (const auto pe_uvw : entityfounds)
         {
            std::unique_ptr<xmapping::xMapping> pmapping(xMappingBuilderHolderSingleton::instance().buildMapping(*pe_uvw.first));
            std::cout << "Entity " << pe_uvw.first << " local Coordinate " << pe_uvw.second << " Global Coordinate "
                      << pmapping->eval(pe_uvw.second) << std::endl;
         }
      }
      std::ofstream out("octree_3d_from_xmesh_getOctreeGrid.pos");
      exportGMSH(mesh.getOctreeGrid(), out);
   }
   return xMPIEnv::finalize();
}
