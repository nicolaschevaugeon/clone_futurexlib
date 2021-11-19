/*
   xfem : C++ Finite Element Library
   developed under the GNU Lesser General Public License
   See the NOTICE & LICENSE files for conditions.
*/

#ifndef _SURF_2_LSET_H_
#define _SURF_2_LSET_H_

#include <limits>

// xoctree
#include "oExport.h"
#include "oKeyManager.h"
#include "oOctree.h"

// xfem
#include "xLevelSet.h"
#include "xNearestNeighborInterface.h"

/// \brief Interface for creating a levelSet from a surface mesh
/// \param[in] mesh_surf : Surface mesh
/// \param[in] bb_ratio : How much the bounding box of the surface mesh is enlarged to create the background mesh. double[3]
/// \param[in] elemPerEdge : For creating a classical regular AOMD mesh (NOT octree): number of elements per edges.
/// \param[in] levelMax : For creating an octree mesh
namespace xinterface
{
namespace xoctree
{
class surface2LevelSet
{
  public:
   surface2LevelSet(xfem::xMesh &mesh_surf_, double *bb_ratio_, double elemSize_, int elemPerEdge_ = -1);
   surface2LevelSet(xfem::xMesh &mesh_surf_, double *bb_ratio_, int levelMax_);
   void createGeoFile(string filename);
   void createMeshBBox(xfem::xMesh *meshbox);
   /// Compute the Lset, export it in a file, and fill a xLevelSet field if given
   void computeLsOnMeshAndExport(xfem::xMesh &mesh, string lsname, xfem::xLevelSet *ls = nullptr);
   /// Compute the Lset, on the octree (finest level), export it to a .ls3d file and fill the lsvVec if given
   void computeLsOnOctreeAndExport(string lsname, std::vector<double> *lsvVec = nullptr);
   void setOctreeLevelMax(int lmax) { levelMax = lmax; }
   const xgeom::xBoundingBox &getBB() const { return BB; }
   void enlargeBB(double scalex, double scaley, double scalez)
   {
      auto center = BB.center();
      auto length = BB.diag();
      length(0) *= scalex;
      length(1) *= scaley;
      length(2) *= scalez;
      BB.min = center - length * 0.5;
      BB.max = center + length * 0.5;
   }

  private:
   xtensor::xVector<> normaleToTri(const AOMD::mEntity *tri);
   xtensor::xVector<> normaleToTri(AOMD::mEntity *vv1, AOMD::mEntity *vv2, AOMD::mEntity *vv3);

   xfem::xMesh &mesh_surf;
   double *bb_ratio, elemSize;
   int levelMax;
   xgeom::xBoundingBox BB;
   xNearestNeighborInterface<xfem::xIter> annint;
};

class computeWeightedPseudoNormals
{
  public:
   computeWeightedPseudoNormals(xfem::xMesh &mesh_) : mesh(mesh_), pi(4. * atan(1.)) {}

   void proceed();

   xtensor::xVector<> getNormale(const AOMD::mEntity *e) const
   {
      // return pseudoNormals[e];
      std::map<const AOMD::mEntity *, xtensor::xVector<>>::const_iterator itf = pseudoNormals.find(e);
      if (itf != pseudoNormals.end())
         return itf->second;
      else
         throw;
      //         return pseudoNormals.find(e)->second;
   }

   xtensor::xVector<> normaleToTri(const AOMD::mEntity *tri)
   {
      //    tri->print();
      AOMD::mVertex *v1 = (AOMD::mVertex *)tri->get(0, 0);
      AOMD::mVertex *v2 = (AOMD::mVertex *)tri->get(0, 1);
      AOMD::mVertex *v3 = (AOMD::mVertex *)tri->get(0, 2);

      xtensor::xVector<> v12(v1->point(), v2->point());
      xtensor::xVector<> v13(v1->point(), v3->point());

      xtensor::xVector<> normale = v12 % v13;
      normale.norm();

      return normale;
   }

  private:
   xfem::xMesh &mesh;
   std::map<const AOMD::mEntity *, xtensor::xVector<>> pseudoNormals;
   const double pi;

   std::pair<AOMD::mVertex *, AOMD::mVertex *> getOtherVertice(AOMD::mEntity *face, AOMD::mVertex *vIn)
   {
      AOMD::mVertex *v0 = (AOMD::mVertex *)face->get(0, 0);
      AOMD::mVertex *v1 = (AOMD::mVertex *)face->get(0, 1);
      AOMD::mVertex *v2 = (AOMD::mVertex *)face->get(0, 2);

      if (vIn == v0)
      {
         return std::make_pair(v1, v2);
      }
      else if (vIn == v1)
      {
         return std::make_pair(v0, v2);
      }
      else
      {
         return std::make_pair(v1, v0);
      }
   }
};

#ifdef HAVE_CGAL
// class surface2LevelSetGeom{

// public:

//  surface2LevelSetGeom(xfem::xMesh &mesh_surf_, double *bb_ratio_, int levelMax_, int recurMax = 0)
//    : mesh_surf(mesh_surf_), bb_ratio(bb_ratio_), elemSize(-1.), levelMax(levelMax_), weightedNormalsComp(mesh_surf_)
//  {

//    cout<<"inside\n";

//    mesh_surf.compute_bounding_box(BBmin,BBmax);
//    cout<<"construc2\n";
//    enlargeBB(bb_ratio_[0],bb_ratio_[1],bb_ratio_[2]);

//    cout<<"compute weighted pseudo-normals...";
//    weightedNormalsComp.proceed();
//    cout<<" done !\n";

////    AOMD::AOMD_Util::Instance()->ex_port("out.msh", &mesh_surf);
//  };

//  /// Compute the Lset, on the octree (finest level), export it to a .ls3d file and fill the lsvVec if given
//  void computeLsOnOctreeAndExport(string lsname, vector<double> *lsvVec=0);
//  void setOctreeLevelMax(int lmax){levelMax = lmax;}
//  void getBB(xtensor::xPoint &BBmin_, xtensor::xPoint &BBmax_) {BBmin_=BBmin; BBmax_=BBmax; return;}
//  void enlargeBB(double scalex, double scaley, double scalez){
//    xtensor::xPoint center=(BBmin+ BBmax)*0.5;
//    xtensor::xPoint length = (BBmax - BBmin);
//    length(0)*=scalex; length(1)*=scaley; length(2)*=scalez;
//    BBmin = center-length*0.5;
//    BBmax = center+length*0.5;
//  }

// private:

//  xtensor::xVector normaleToTri(const AOMD::mEntity *tri);
//  xtensor::xVector normaleToTri(AOMD::mEntity *vv1, AOMD::mEntity *vv2, AOMD::mEntity *vv3);
//  std::pair<int, int> getSubEntity(AOMD::mEntity *nearestEntity, xtensor::xPoint &closestPointOnsurf, bool verbose=false);

// private:
//  xfem::xMesh &mesh_surf;
//  double *bb_ratio, elemSize;
//  int levelMax;
//  xtensor::xPoint BBmin, BBmax;
////  std::map<AOMD::mVertex *, xtensor::xVector> mapNormales;
//  computeWeightedPseudoNormals weightedNormalsComp;
//};

class surface2LevelSetNarrowBand
{
  public:
   surface2LevelSetNarrowBand(xfem::xMesh &mesh_surf_, double *bb_ratio_, int levelMax_, bool useWeightedPseudoNormals_ = true,
                              bool cubicScale = false);
   /// Compute lset on narrow band
   void computeLsOnNarrowBand(std::vector<double> *lsvVec);
   /// Extend the levelset on the whole domain using fastmarching like algorithm
   void extendLset(std::vector<double> *lsvVec);
   //// Save levelset to file (ascii)
   void saveLset(std::string lsname, std::vector<double> *lsvVec);

   /// returns bounding box
   void getBB(xtensor::xPoint &BBmin_, xtensor::xPoint &BBmax_)
   {
      BBmin_ = BBmin;
      BBmax_ = BBmax;
      return;
   }
   /// enlarge bouding box from x, y and z scales
   void enlargeBB(double scalex, double scaley, double scalez)
   {
      xtensor::xPoint center = (BBmin + BBmax) * 0.5;
      xtensor::xPoint length = (BBmax - BBmin);
      length(0) *= scalex;
      length(1) *= scaley;
      length(2) *= scalez;
      BBmin = center - length * 0.5;
      BBmax = center + length * 0.5;
   }

   ~surface2LevelSetNarrowBand()
   {
      if (poctree) delete poctree;
      if (pmapping) delete pmapping;
   }

   /// Activate cells close to the surface mesh : iEpsilon = nb de cellules ajoutees en plus de chaque cote (marge...)
   void activateBoundingCells(int iEpsilon = 1);
   void activateAllCells();

   void releaseOctree()
   {
      if (poctree) delete poctree;
      if (pmapping) delete pmapping;
   }

  private:
   /// Gives wheter a point is on the face, one of its edge (and which one), or on one of its vertices (which one ?)
   std::pair<int, int> getSubEntity(AOMD::mEntity *nearestEntity, const xtensor::xPoint &closestPointOnsurf,
                                    bool verbose = false);

   /// Change the bounding box to make it cubic
   void changeRatioToCubic(double *ratios);

   xfem::xMesh &mesh_surf;
   double *bb_ratio, elemSize;
   int levelMax;
   xtensor::xPoint BBmin, BBmax;
   std::map<AOMD::mVertex *, xtensor::xVector<>> mapNormales;
   computeWeightedPseudoNormals weightedNormalsComp;
   ::xoctree::oTopo topo;
   ::xoctree::oMappingCartesian *pmapping;
   ::xoctree::oOctree *poctree;
   bool useWeightedPseudoNormals;
};
#endif

}  // namespace xoctree
}  // namespace xinterface

#endif
