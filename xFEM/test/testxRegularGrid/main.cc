#include <fstream>

#include "xMPIEnv.h"
#include "xMappingBuilderHolder.h"
#include "xMesh.h"
#include "xOctreeGrid.h"

using namespace xfem;
using namespace xtensor;
using namespace xtool;

int main(int argc, char* argv[])
{
   xMPIEnv::init(argc, argv);

   {  // A case on a 3d mesh, using the xMesh internal createOctree and locateOctree
      xMesh mesh("data/CubeWithCubeHole.msh");
      std::cout << "Mesh read : " << mesh.size(0) << " vertices " << mesh.size(3) << " solids " << std::endl;
      std::array<xPoint, 2> targets{xPoint{0.05, 0.05, 0.05}, xPoint{0., 0., 0.}};
      std::cout << "Testing mesh.locateElement" << std::endl;
      for (auto target : targets)
      {
         auto entityfounds = mesh.locateElement(target);
         std::cout << "Point " << target << " Found in " << std::distance(entityfounds.begin(), entityfounds.end())
                   << " Elements " << std::endl;
         for (const auto pe_uvw : entityfounds)
         {
            std::unique_ptr<xmapping::xMapping> pmapping(xMappingBuilderHolderSingleton::instance().buildMapping(*pe_uvw.first));
            std::cout << "Entity " << pe_uvw.first << " local Coordinate " << pe_uvw.second << " Global Coordinate "
                      << pmapping->eval(pe_uvw.second) << std::endl;
         }
      }
      std::ofstream out("regular_grid_from_xmesh_getRegularGrid.txt");
      mesh.getRegularGrid().print(out);
   }
   return xMPIEnv::finalize();
}
