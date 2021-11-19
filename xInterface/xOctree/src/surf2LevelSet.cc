/*
   xfem : C++ Finite Element Library
   developed under the GNU Lesser General Public License
   See the NOTICE & LICENSE files for conditions.
*/

#include "surf2LevelSet.h"

#include <fstream>
#include <iostream>

#include "oOctree.h"

#ifdef HAVE_CGAL
// Geom
#include "xDistanceNearestPoint.h"
#include "xDistanceNearestPointGenerator.h"
#endif

#include "xMesh.h"

using AOMD::mEntity;
using AOMD::mVertex;

using xfem::xMesh;
using xgeom::xBoundingBox;
using ::xoctree::ExportGMSHAsciiNodes;
using ::xoctree::oKey;
using ::xoctree::oKeyManager;
using ::xoctree::oMappingCartesian;
using ::xoctree::oOctree;
using xtensor::xPoint;
using namespace std;
using namespace xfem;

namespace xinterface
{
namespace xoctree
{
surface2LevelSet::surface2LevelSet(xMesh &mesh_surf_, double *bb_ratio_, double elemSize_, int elemPerEdge_)
    : mesh_surf(mesh_surf_), bb_ratio(bb_ratio_), elemSize(elemSize_), levelMax(-1), annint(mesh_surf.begin(0), mesh_surf.end(0))
{
   BB = mesh_surf.compute_bounding_box();
   enlargeBB(bb_ratio[0], bb_ratio[1], bb_ratio[2]);

   if (elemPerEdge_ != -1)
   {
      auto lenght = BB.diag();
      double smallerSize = min(min(lenght[0], lenght[1]), lenght[2]);
      elemSize = smallerSize / elemPerEdge_;
   }
};

surface2LevelSet::surface2LevelSet(xMesh &mesh_surf_, double *bb_ratio_, int levelMax_)
    : mesh_surf(mesh_surf_), bb_ratio(bb_ratio_), elemSize(-1.), levelMax(levelMax_), annint(mesh_surf.begin(0), mesh_surf.end(0))
{
   BB = mesh_surf.compute_bounding_box();
   enlargeBB(bb_ratio_[0], bb_ratio_[1], bb_ratio_[2]);
};

xtensor::xVector<> surface2LevelSet::normaleToTri(const AOMD::mEntity *tri)
{
   // tri->print();
   mVertex *v1 = (mVertex *)tri->get(0, 0);
   mVertex *v2 = (mVertex *)tri->get(0, 1);
   mVertex *v3 = (mVertex *)tri->get(0, 2);

   xtensor::xVector<> v12(v1->point(), v2->point());
   xtensor::xVector<> v13(v1->point(), v3->point());

   xtensor::xVector<> normale = v12 % v13;
   normale.norm();

   //  cout<<"norm="<<normale<<endl;

   return normale;
}

xtensor::xVector<> surface2LevelSet::normaleToTri(mEntity *vv1, mEntity *vv2, mEntity *vv3)
{
   mVertex *v1 = (mVertex *)vv1;
   mVertex *v2 = (mVertex *)vv2;
   mVertex *v3 = (mVertex *)vv3;

   xtensor::xVector<> v12(v1->point(), v2->point());
   xtensor::xVector<> v13(v1->point(), v3->point());

   xtensor::xVector<> normale = v12 % v13;
   normale.norm();

   return normale;
}

void surface2LevelSet::createGeoFile(string filename)
{
   fstream file(filename.c_str(), ios::in | ios::out | ios::trunc);

   file << "x0 = " << BB.min(0) << ";" << endl;
   file << "y0 = " << BB.min(1) << ";" << endl;
   file << "z0 = " << BB.min(2) << ";" << endl;

   file << "x1 = " << BB.max(0) << ";" << endl;
   file << "y1 = " << BB.max(1) << ";" << endl;
   file << "z1 = " << BB.max(2) << ";" << endl;

   file << "h = " << elemSize << ";" << endl;

   file << "Point(1) = {x0,y0,z0,h};" << endl;
   file << "Point(2) = {x1,y0,z0,h};" << endl;
   file << "Point(3) = {x1,y1,z0,h};" << endl;
   file << "Point(4) = {x0,y1,z0,h};" << endl;
   file << "Point(5) = {x0,y0,z1,h};" << endl;
   file << "Point(6) = {x1,y0,z1,h};" << endl;
   file << "Point(7) = {x1,y1,z1,h};" << endl;
   file << "Point(8) = {x0,y1,z1,h};" << endl;

   file << "Line(1) = {1,2};" << endl;
   file << "Line(2) = {2,3};" << endl;
   file << "Line(3) = {3,4};" << endl;
   file << "Line(4) = {4,1};" << endl;
   file << "Line(5) = {1,5};" << endl;
   file << "Line(6) = {2,6};" << endl;
   file << "Line(7) = {3,7};" << endl;
   file << "Line(8) = {4,8};" << endl;
   file << "Line(9) = {5,6};" << endl;
   file << "Line(10) = {6,7};" << endl;
   file << "Line(11) = {7,8};" << endl;
   file << "Line(12) = {8,5};" << endl;

   file << "Line Loop(13) = {1,2,3,4};" << endl;
   file << "Plane Surface(20) = {13};" << endl;
   file << "Line Loop(15) = {-9,-5,1,6};" << endl;
   file << "Plane Surface(16) = {15};" << endl;
   file << "Line Loop(17) = {-10,-6,2,7};" << endl;
   file << "Plane Surface(18) = {17};" << endl;
   file << "Line Loop(19) = {-11,-7,3,8};" << endl;
   file << "Plane Surface(22) = {19};" << endl;
   file << "Line Loop(21) = {12,-5,-4,8};" << endl;
   file << "Plane Surface(14) = {21};" << endl;
   file << "Line Loop(23) = {11,12,9,10};" << endl;
   file << "Plane Surface(24) = {23};" << endl;

   file << "Surface Loop(25) = {14,16,18,20,22,24};" << endl;
   file << "Volume(26) = {25};" << endl;

   file << "Physical Point   (101)  = {1} ;" << endl;
   file << "Physical Point   (102)  = {2} ;" << endl;
   file << "Physical Point   (103)  = {3} ;" << endl;

   file << "Physical Surface   (214) = {14};" << endl;
   file << "Physical Surface   (109) = {16};" << endl;
   file << "Physical Surface   (218) = {18};" << endl;
   file << "Physical Surface   (220) = {20};" << endl;
   file << "Physical Surface   (111) = {22};" << endl;
   file << "Physical Surface   (224) = {24};" << endl;

   file << "Physical Volume    (121) = {26};" << endl;

   file.close();
}

inline int linearEncoding(int i, int j, int k, int imax, int jmax, int kmax)
{
   return i + (imax + 1) * j + (imax + 1) * (jmax + 1) * k;
}

void surface2LevelSet::createMeshBBox(xMesh *meshbox)
{
   AOMD::mMesh *mmeshbox = &meshbox->getMesh();
   cout << elemSize << endl;
   auto lenght = BB.diag();
   int elemPerEdgeX = lenght[0] / elemSize;
   int elemPerEdgeY = lenght[1] / elemSize;
   int elemPerEdgeZ = lenght[2] / elemSize;

   double elemSizeX = lenght[0] / elemPerEdgeX;
   double elemSizeY = lenght[1] / elemPerEdgeY;
   double elemSizeZ = lenght[2] / elemPerEdgeZ;

   cout << elemPerEdgeX << " " << elemPerEdgeY << " " << elemPerEdgeZ << endl;

   int id = 0;

   for (int k = 0; k < elemPerEdgeZ + 1; ++k)
   {
      for (int i = 0; i < elemPerEdgeX + 1; ++i)
      {
         for (int j = 0; j < elemPerEdgeY + 1; ++j)
         {
            double xyz[3] = {BB.min(0) + i * elemSizeX, BB.min(1) + j * elemSizeY, BB.min(2) + k * elemSizeZ};
            mmeshbox->createVertex(id, xyz[0], xyz[1], xyz[2], nullptr);
            ++id;
         }
      }
   }

   cout << "Id is" << id << endl;
   const int dim = 3;
   for (int i = 0; i < elemPerEdgeX; ++i)
   {
      for (int j = 0; j < elemPerEdgeY; ++j)
      {
         for (int k = 0; k < elemPerEdgeZ; ++k)
         {
            int ids[8] = {linearEncoding(i, j, k, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ),
                          linearEncoding(i + 1, j, k, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ),
                          linearEncoding(i + 1, j + 1, k, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ),
                          linearEncoding(i, j + 1, k, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ),
                          linearEncoding(i, j, k + 1, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ),
                          linearEncoding(i + 1, j, k + 1, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ),
                          linearEncoding(i + 1, j + 1, k + 1, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ),
                          linearEncoding(i, j + 1, k + 1, elemPerEdgeX, elemPerEdgeY, elemPerEdgeZ)};

            //        cout<<"Ids=";
            //        for(int kk=0;kk<8;++kk){
            //          cout<<ids[kk]<<" ";
            //        }
            //        cout<<endl;

#if 1
            mmeshbox->createHexWithVertices(ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7],
                                            mmeshbox->getGEntity(100, dim));
#else
            // Maillage non alterne
            mEntity *tet = 0;
            tet = (mEntity *)mesh.createTetWithVertices(ids[0], ids[1], ids[5], ids[3], mesh.getGEntity(100, 3));
            tet = (mEntity *)mesh.createTetWithVertices(ids[1], ids[2], ids[3], ids[5], mesh.getGEntity(100, 3));
            tet = (mEntity *)mesh.createTetWithVertices(ids[2], ids[3], ids[5], ids[6], mesh.getGEntity(100, 3));
            tet = (mEntity *)mesh.createTetWithVertices(ids[0], ids[3], ids[4], ids[5], mesh.getGEntity(100, 3));
            tet = (mEntity *)mesh.createTetWithVertices(ids[3], ids[4], ids[5], ids[7], mesh.getGEntity(100, 3));
            tet = (mEntity *)mesh.createTetWithVertices(ids[3], ids[5], ids[6], ids[7], mesh.getGEntity(100, 3));
#endif
         }
      }
   }

   classifyUnclassifiedVerices(mmeshbox);
   // meshbox->modifyAllState();
   mmeshbox->modifyState(3, 2, true);
   mmeshbox->modifyState(3, 1, true);
   mmeshbox->modifyState(3, 0, true);
   mmeshbox->modifyState(2, 1, true);
   mmeshbox->modifyState(2, 0, true);
   mmeshbox->modifyState(1, 0, true);
   mmeshbox->modifyState(0, 1, true);
   mmeshbox->modifyState(0, 2, true);
   mmeshbox->modifyState(0, 3, true);
   mmeshbox->modifyState(1, 2, true);
   mmeshbox->modifyState(1, 3, true);
   mmeshbox->modifyState(2, 3, true);
   AOMD::AOMD_Util::Instance()->ex_port("meshBox.msh", mmeshbox);

   //  Export au format mshV2:
   // TODO!!!
}

void surface2LevelSet::computeLsOnMeshAndExport(xMesh &mesh, string lsname, xLevelSet *ls)
{
   fstream file(lsname.c_str(), ios::in | ios::out | ios::trunc);

   int dim = 3;

   file << "$NodeData" << endl;
   file << "1" << endl << "\"" << lsname << "\"" << endl;
   file << "1" << endl << "0.0" << endl;
   file << "3" << endl << "0" << endl << "1" << endl;
   file << mesh.size(0) << endl;

   for (xIter it = mesh.begin(0); it != mesh.end(0); ++it)
   {
      mVertex *v = (mVertex *)*it;
      double d = 0.;
      mVertex *closestonsurf = annint.nearestvertexto(*v, d);

      xtensor::xVector<> normale, normLoc;

      for (int iv = 0; iv < closestonsurf->size(dim - 1); ++iv)
      {
         mEntity *voisin = closestonsurf->get(dim - 1, iv);
         xtensor::xVector<> normLoc = normaleToTri(voisin);
         normale += normLoc;
      }

      normale.norm();

      xtensor::xVector<> direction(closestonsurf->point(), v->point());
      double dir = direction * normale;
      if (dir < 0.)
         dir = -1.;
      else
         dir = 1.;

      file << v->getId() << " " << dir * d << endl;
      if (ls) ls->operator()(v) = dir * d;
   }

   file << "$EndNodeData" << endl;

   file.close();
}

void surface2LevelSet::computeLsOnOctreeAndExport(string lsname, vector<double> *lsvVec)
{
   if (levelMax == -1) throw;

   fstream file(lsname.c_str(), ios::in | ios::out | ios::trunc);

   int dim = 3;
   // int sizeMax = (pow(2.,levelMax)+1)* (pow(2.,levelMax)+1) * (pow(2.,levelMax)+1);

   int imax = pow(2., levelMax) + 1;
   int jmax = pow(2., levelMax) + 1;
   int kmax = pow(2., levelMax) + 1;

   if (lsvVec) lsvVec->reserve(imax * jmax * kmax);
   const auto lenght = BB.diag();
   double hx = lenght[0] / (imax - 1);
   double hy = lenght[1] / (jmax - 1);
   double hz = lenght[2] / (kmax - 1);

   file << dim << endl;  // dim
   //  file<<(BBmax(0)-BBmin(0))<<endl;//sizeX
   //  file<<(BBmax(1)-BBmin(1))<<endl;//sizeY
   //  file<<(BBmax(2)-BBmin(2))<<endl;//sizeZ
   file << imax << endl;       // sizeX
   file << jmax << endl;       // sizeY
   file << kmax << endl;       // sizeZ
   file << hx << endl;         // spaceX
   file << hy << endl;         // spaceY
   file << hz << endl;         // spaceZ
   file << BB.min(0) << endl;  // originX
   file << BB.min(1) << endl;  // originY
   file << BB.min(2) << endl;  // originZ

   for (int k = 0; k < kmax; ++k)
   {
      for (int j = 0; j < jmax; ++j)
      {
         for (int i = 0; i < imax; ++i)
         {  // 1 noeud de plus que de cells

            // int ijk[3]={i,j,k};
            double P[3] = {BB.min(0) + i * hx, BB.min(1) + j * hy, BB.min(2) + k * hz};

            double d;
            mVertex *closestonsurf = annint.nearestvertexto(P, d);

            xtensor::xVector<> direction(closestonsurf->point(), xPoint(P[0], P[1], P[2]));
            xtensor::xVector<> normale(0., 0., 0.);

            for (int iv = 0; iv < closestonsurf->size(dim - 1); ++iv)
            {
               mEntity *voisin = closestonsurf->get(dim - 1, iv);
               xtensor::xVector<> normLoc = normaleToTri(voisin);
               normale += normLoc;
            }

            normale.norm();

            double dir = direction * normale;
            if (dir < 0.)
               dir = -1.;
            else
               dir = 1.;

            file << dir * d << endl;
            if (lsvVec) lsvVec->push_back(dir * d);
         }
      }
   }

   file.close();
}

void computeWeightedPseudoNormals::proceed()
{
   // Faces
   for (xIter it = mesh.begin(2); it != mesh.end(2); ++it)
   {
      AOMD::mEntity *tri = *it;
      pseudoNormals[tri] = 2. * pi * normaleToTri(tri);
   }

   // Edges
   for (xIter it = mesh.begin(1); it != mesh.end(1); ++it)
   {
      AOMD::mEntity *edge = *it;

      // Valable seulement en 3D
      if (mesh.size(2))
      {
         // Pour traiter le cas ou la surface n'est pas fermee : on peut avoir des edges avec une seule facette
         AOMD::mEntity *face1 = edge->get(2, 0);

         if (edge->size(2) > 1)
         {
            AOMD::mEntity *face2 = edge->get(2, 1);
            pseudoNormals[edge] = pi * (normaleToTri(face1) + normaleToTri(face2));
         }
         else
         {
            pseudoNormals[edge] = normaleToTri(face1);
         }
      }
      else
      {
         //            cout<<"WeightedNormals in 1D\n";
         // En 2D:
         pseudoNormals[edge] = (xtensor::xVector<>(static_cast<mVertex *>(edge->get(0, 0))->point(),
                                                   static_cast<mVertex *>(edge->get(0, 1))->point()) %
                                xtensor::xVector<>(0, 0, 1))
                                   .norm();

         //            edge->print();
         //            cout<<pseudoNormals[edge]<<endl;
      }
   }

   // Vertex
   for (xIter it = mesh.begin(0); it != mesh.end(0); ++it)
   {
      AOMD::mVertex *v = (AOMD::mVertex *)*it;
      xtensor::xVector<> normale(0, 0, 0);

      for (int iface = 0; iface < v->size(2); ++iface)
      {
         AOMD::mEntity *face = v->get(2, iface);

         std::pair<AOMD::mVertex *, AOMD::mVertex *> vOthers = getOtherVertice(face, v);

         xtensor::xVector<> vA(v->point(), vOthers.first->point());
         xtensor::xVector<> vB(v->point(), vOthers.second->point());
         vA.norm();
         vB.norm();
         double angle = acos(vA * vB);

         normale += normaleToTri(face) * angle;
      }

      pseudoNormals[v] = normale;
   }
}

#ifdef HAVE_CGAL
//--------------------------------------
surface2LevelSetNarrowBand::surface2LevelSetNarrowBand(xfem::xMesh &mesh_surf_, double *bb_ratio_, int levelMax_,
                                                       bool useWeightedPseudoNormals_, bool cubicScale)
    : mesh_surf(mesh_surf_),
      bb_ratio(bb_ratio_),
      elemSize(-1.),
      levelMax(levelMax_),
      weightedNormalsComp(mesh_surf_),
      pmapping(nullptr),
      poctree(nullptr),
      useWeightedPseudoNormals(useWeightedPseudoNormals_)
{
   cout << "inside\n";

   mesh_surf.compute_bounding_box(BBmin, BBmax);
   cout << "construc2\n";

   if (cubicScale) changeRatioToCubic(bb_ratio);
   enlargeBB(bb_ratio[0], bb_ratio[1], bb_ratio[2]);

   cout << "Compute weighted pseudo-normals on surface mesh...";
   if (useWeightedPseudoNormals) weightedNormalsComp.proceed();
   cout << " done !\n";

   cout << "Create local octree...";
   double box_inf[3] = {BBmin(0), BBmin(1), BBmin(2)};
   double box_sup[3] = {BBmax(0), BBmax(1), BBmax(2)};
   int per[3] = {0, 0, 0};
   pmapping = new oMappingCartesian(topo, 3, box_inf, box_sup);
   poctree = new oOctree(*pmapping, levelMax, per);
   poctree->deactivateFinestLevel();
   cout << " done !\n";
}

#if 0

    xtensor::xVector surface2LevelSetGeom::normaleToTri(const mEntity *tri){
      //tri->print();
      mVertex *v1= (mVertex *) tri->get(0,0);
      mVertex *v2= (mVertex *) tri->get(0,1);
      mVertex *v3= (mVertex *) tri->get(0,2);

      xtensor::xVector v12(v1->point(),v2->point());
      xtensor::xVector v13(v1->point(),v3->point());

      xtensor::xVector normale = v12%v13;
      normale.norm();

      return normale;

    }

    xtensor::xVector surface2LevelSetGeom::normaleToTri(mEntity *vv1, mEntity *vv2, mEntity *vv3){

      mVertex *v1= (mVertex *) vv1;
      mVertex *v2= (mVertex *) vv2;
      mVertex *v3= (mVertex *) vv3;

      xtensor::xVector v12(v1->point(),v2->point());
      xtensor::xVector v13(v1->point(),v3->point());

      xtensor::xVector normale = v12%v13;
      normale.norm();

      return normale;

    }




    void surface2LevelSetGeom::computeLsOnMeshAndExport(xMesh &mesh, string lsname, xLevelSet *ls){

      //  fstream file(lsname.c_str(), ios::in | ios::out | ios::trunc);

      //  int dim=3;

      //  file<<"$NodeData"<<endl;
      //  file<<"1"<<endl<<"\""<<lsname<<"\""<<endl;
      //  file<<"1"<<endl<<"0.0"<<endl;
      //  file<<"3"<<endl<<"0"<<endl<<"1"<<endl;
      //  file<<mesh.size(0)<<endl;




      //  for(xIter it=mesh.begin(0); it!=mesh.end(0); ++it){

      //    mVertex *v= (mVertex *) *it;
      //    double d=0.;
      //    mVertex *closestonsurf = annint.nearestvertexto(*v, d );

      //    xtensor::xVector normale, normLoc;

      //    for(int iv=0; iv<closestonsurf->size(dim-1); ++iv){
      //      mEntity *voisin=closestonsurf->get(dim-1,iv);
      //      xtensor::xVector normLoc=normaleToTri(voisin);
      //      normale+=normLoc;
      //    }

      //    normale.norm();

      //    xtensor::xVector direction(closestonsurf->point(),v->point());
      //    double dir= direction*normale;
      //    if(dir<0.) dir=-1.;
      //    else dir=1.;

      //    file<<v->getId()<<" "<< dir*d<<endl;
      //    if(ls) ls->operator()(v)=dir*d;

      //  }

      //  file<<"$EndNodeData"<<endl;

      //  file.close();
    }

    void surface2LevelSetGeom::computeLsOnOctreeAndExport(string lsname, vector<double> *lsvVec){

      //    xNearestNeighborInterface<xIter > annint(mesh_surf.begin(0), mesh_surf.end(0));
      //    xgeom::xDistanceNearestPoint < xIter, xgeom::xAOMDToANN > distance_nearest_point(mesh_surf.begin(0), mesh_surf.end(0));
      xgeom::xDistanceNearestPoint < xIter, xgeom::xAOMDTriToCGALTri > distance_nearest_point(mesh_surf.begin(2), mesh_surf.end(2));
      //    xgeom::xDistanceNearestPoint < xIter, xgeom::xAOMDEdgeToCGALEdge > distance_nearest_point(mesh_surf.begin(1), mesh_surf.end(1));
      //    xgeom::xDistanceNearestPoint < xIter, xgeom::xAOMDVertexToCGALPoint > distance_nearest_point(mesh_surf.begin(0), mesh_surf.end(0));


      if(levelMax==-1) throw;

      fstream file(lsname.c_str(), ios::in | ios::out | ios::trunc);

      int dim=3;
      int sizeMax = (int) ((pow(2.,levelMax)+1)* (pow(2.,levelMax)+1) * (pow(2.,levelMax)+1));

      int imax = (int) pow(2.,levelMax)+1;
      int jmax = (int) pow(2.,levelMax)+1;
      int kmax = (int) pow(2.,levelMax)+1;

      if(lsvVec) lsvVec->reserve(imax*jmax*kmax);

      double hx=(BBmax(0)-BBmin(0))/(imax-1);
      double hy=(BBmax(1)-BBmin(1))/(jmax-1);
      double hz=(BBmax(2)-BBmin(2))/(kmax-1);


      file<<dim<<endl;//dim
      //  file<<(BBmax(0)-BBmin(0))<<endl;//sizeX
      //  file<<(BBmax(1)-BBmin(1))<<endl;//sizeY
      //  file<<(BBmax(2)-BBmin(2))<<endl;//sizeZ
      file<<imax<<endl;//sizeX
      file<<jmax<<endl;//sizeY
      file<<kmax<<endl;//sizeZ
      file<<hx<<endl;//spaceX
      file<<hy<<endl;//spaceY
      file<<hz<<endl;//spaceZ
      file<<BBmin(0)<<endl;//originX
      file<<BBmin(1)<<endl;//originY
      file<<BBmin(2)<<endl;//originZ


      double d;
      mEntity *nearestEntity;
      xPoint closestPointOnsurf;

      bool verbose(0);

      for(int k=0; k< kmax; ++k){
	for(int j=0; j< jmax; ++j){
          for(int i=0; i< imax; ++i){//1 noeud de plus que de cells

	    //afficher ijk
	    //          cout<<"ijk="<<i<<" "<<j<<" "<<k<<" "<<i+33*j+33*33*k<<endl;
	    //          if(i+33*j+33*33*k == 4422) verbose=1;
	    if(i+imax*j+imax*jmax*k == 1998) verbose=1;
	    if(verbose) cout<<"id node="<<i+33*j+33*33*k<<" ijk="<<i<<" "<<j<<" "<<k<<endl;

	    int ijk[3]={i,j,k};
	    double P[3]={BBmin(0)+i*hx,BBmin(1)+j*hy,BBmin(2)+k*hz};

	    //        mVertex *closestonsurf = annint.nearestvertexto(P, d );
        mVertex V(0,mPoint(P[0],P[1],P[2]),0);
        if(verbose) cout<<"gridpoint "<<xPoint(P[0],P[1],P[2])<<endl;

	    distance_nearest_point.nearestEntityDistance(&V,nearestEntity,closestPointOnsurf,d);
	    //        d=distance_nearest_point.distance(&V);
	    //        distance_nearest_point.nearestPointDistance(&V, closestPointOnsurf, d);

	    if(verbose) cout<<"nearestpoint "<<closestPointOnsurf<<endl;

#if 1
	    xtensor::xVector direction(closestPointOnsurf,V.point());
	    //        xtensor::xVector normale = normaleToTri(nearestEntity);
	    std::pair<int, int> subEntity = getSubEntity(nearestEntity, closestPointOnsurf, verbose);
	    if(verbose) cout<<subEntity.first<<" "<<subEntity.second<<endl;
	    if(verbose) nearestEntity->print();

	    xtensor::xVector normale;
	    if(subEntity.first <2) normale = weightedNormalsComp.getNormale(nearestEntity->get(subEntity.first, subEntity.second));
	    else normale = weightedNormalsComp.getNormale(nearestEntity);

	    if(normale.mag() < 1.e-12) verbose=1;

	    if(verbose) cout<<"norm "<<normale<<" direction "<< direction <<endl;


	    if(normale.mag() < 1.e-12){
	      //Cas pathologique de normale quasi nulle...
	      //On prend la normale du segement ?
	      if(verbose) cout<<"vanishing normal..\n";
	      if(subEntity.first == 1) throw;//Not coded for vanishing edge normals

	      normale = weightedNormalsComp.getNormale(nearestEntity->get(subEntity.first+1, subEntity.second));

	      if(verbose) cout<<"New norm "<<normale<<endl;

	    }

	    //            normale = mapNormales[closestonsurf];


	    double dir= direction*normale;
	    if(verbose) cout<<"dir "<<dir<<endl;
	    if(dir<0.) dir=-1.;
	    else dir=1.;

	    file<<dir*d<<endl;
	    if(lsvVec) lsvVec->push_back(dir*d);
	    verbose=0;
#endif
	  }
        }
      }

      file.close();
    }


    std::pair<int, int> surface2LevelSetGeom::getSubEntity(mEntity *nearestEntity, xPoint &closestPointOnsurf, bool verbose){

      xfem::xGeomElem geo(nearestEntity);
      geo.setUVWForXYZ(closestPointOnsurf);
      xPoint uvw = geo.getUVW();

      double u(abs(uvw(0))), v(abs(uvw(1)));
      //double epsilon = 2. * std::numeric_limits<double>::epsilon();
      double epsilon = 1.e-13;
      if(verbose) cout<<"epsilon="<<epsilon<<endl;
      if(verbose) cout<<"u="<<u<<" v="<<v<<endl;

      if(u >=epsilon && v >= epsilon && u+v <= 1. -epsilon){
	return make_pair(2,0);
      }else if(u < epsilon){
	if(v < epsilon) return make_pair(0,0);
	if(v > 1. - epsilon) return make_pair(0,2);
	return make_pair(1,2);
      }else if(v < epsilon){
	if(u < epsilon) return make_pair(0,0);
	if(u > 1. - epsilon) return make_pair(0,1);
	return make_pair(1,0);
      }else if( u+v > 1. - epsilon) {
	if(u < epsilon) return make_pair(0,2);
	if(u > 1. - epsilon) return make_pair(0,1);
	return make_pair(1,1);
      }

      cout<<"HaraKiri\n";
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void surface2LevelSetGeom::fillNodesForDistance(int i, int recurMax, xPoint p1, xPoint p2, xPoint p3, std::vector<xPoint> &containeurPts ){
      if(i == recurMax) return ;
      else {

	//Calculer les nouveaux points et les ajouter au containeur
    xPoint p12((p1(0)+p2(0))/2,(p1(1)+p2(1))/2,(p1(2)+p2(2))/2);
    xPoint p13((p1(0)+p3(0))/2,(p1(1)+p3(1))/2,(p1(2)+p3(2))/2);
    xPoint p23((p3(0)+p2(0))/2,(p3(1)+p2(1))/2,(p3(2)+p2(2))/2);

    xPoint newP1((p1(0)+p12(0)+p13(0))/3,(p1(1)+p12(1)+p13(1))/3,(p1(2)+p12(2)+p13(2))/3);
	containeurPts.push_back(newP1);

    xPoint newP2((p2(0)+p12(0)+p23(0))/3,(p2(1)+p12(1)+p23(1))/3,(p2(2)+p12(2)+p23(2))/3);
	containeurPts.push_back(newP2);

    xPoint newP3((p3(0)+p13(0)+p23(0))/3,(p3(1)+p13(1)+p23(1))/3,(p3(2)+p13(2)+p23(2))/3);
	containeurPts.push_back(newP3);

	if (i==0){
          xPoint newP4((p12(0)+p13(0)+p23(0))/3,(p12(1)+p13(1)+p23(1))/3,(p12(2)+p13(2)+p23(2))/3);
          containeurPts.push_back(newP4);
        }

	//recursivite
	this->fillNodesForDistance(i+1,recurMax,p1,p12,p13,containeurPts);
	this->fillNodesForDistance(i+1,recurMax,p12,p2,p23,containeurPts);
	this->fillNodesForDistance(i+1,recurMax,p3,p13,p23,containeurPts);
	this->fillNodesForDistance(i+1,recurMax,p12,p23,p13,containeurPts);

      }



      return;

    }
#endif

//--------------------------------------
#if 1

void surface2LevelSetNarrowBand::changeRatioToCubic(double *ratios)
{
   cout << "Update BBox dimension ratios to cubic\n";

   // Length of the sides of the BBox along x, y and z
   double sideLength[3] = {abs(BBmin(0) - BBmax(0)), abs(BBmin(1) - BBmax(1)), abs(BBmin(2) - BBmax(2))};

   // Classify decreasing length : classDim[0] = k --> bigger length along direction k
   int classDim[3] = {-1, -1, -1};

   if (sideLength[0] > sideLength[1] && sideLength[0] > sideLength[2])
   {
      // bigX
      classDim[0] = 0;
      if (sideLength[1] > sideLength[2])
      {
         classDim[1] = 1;
         classDim[2] = 2;
      }
      else
      {
         classDim[1] = 2;
         classDim[2] = 1;
      }
   }

   if (sideLength[1] > sideLength[0] && sideLength[1] > sideLength[2])
   {
      // bigY
      classDim[0] = 1;
      if (sideLength[0] > sideLength[2])
      {
         classDim[1] = 0;
         classDim[2] = 2;
      }
      else
      {
         classDim[1] = 2;
         classDim[2] = 0;
      }
   }

   if (sideLength[2] > sideLength[1] && sideLength[2] > sideLength[0])
   {
      // bigZ
      classDim[0] = 2;
      if (sideLength[1] > sideLength[0])
      {
         classDim[1] = 1;
         classDim[2] = 0;
      }
      else
      {
         classDim[1] = 0;
         classDim[2] = 1;
      }
   }

   // Scale everything (except the bigger one)
   ratios[classDim[1]] = ratios[classDim[0]] * (sideLength[classDim[0]] / sideLength[classDim[1]]);
   ratios[classDim[2]] = ratios[classDim[0]] * (sideLength[classDim[0]] / sideLength[classDim[2]]);

   printf("New cubic ratios : [%f, %f, %f]", ratios[0], ratios[1], ratios[2]);

   return;
}

void surface2LevelSetNarrowBand::activateBoundingCells(int iEpsilon)
{
   for (mEntity *e : mesh_surf.range(2))
   {
      xgeom::xBoundingBox bbox = xfem::compute_bounding_box(*e);
      xPoint pmin = bbox.min;
      xPoint pmax = bbox.max;

      //        cout<<"bbmin ="<<pmin<<endl;
      //        cout<<"bbmax ="<<pmax<<endl;

      int ijkmin[3] = {0, 0, 0};
      int ijkmax[3] = {0, 0, 0};
      poctree->locatePointOnFinestLevel(pmin(0), pmin(1), pmin(2), ijkmin);
      poctree->locatePointOnFinestLevel(pmax(0), pmax(1), pmax(2), ijkmax);

      //        cout<<"ijkmin = "<<ijkmin[0]<<" "<<ijkmin[1]<<" "<<ijkmin[2]<<endl;
      //        cout<<"ijkmax = "<<ijkmax[0]<<" "<<ijkmax[1]<<" "<<ijkmax[2]<<endl;

      int extent[3] = {ijkmax[0] - ijkmin[0] + 2 * iEpsilon, ijkmax[1] - ijkmin[1] + 2 * iEpsilon,
                       ijkmax[2] - ijkmin[2] + 2 * iEpsilon};

      int imax = (int)pow(2., levelMax) + 1;
      int jmax = (int)pow(2., levelMax) + 1;
      int kmax = (int)pow(2., levelMax) + 1;

      int ijkcurr[3];
      for (int k = 0; k <= extent[2]; ++k)
      {
         for (int j = 0; j <= extent[1]; ++j)
         {
            for (int i = 0; i <= extent[0]; ++i)
            {
               ijkcurr[0] = min(max(ijkmin[0] - iEpsilon + i, 0), imax);
               ijkcurr[1] = min(max(ijkmin[1] - iEpsilon + j, 0), jmax);
               ijkcurr[2] = min(max(ijkmin[2] - iEpsilon + k, 0), kmax);

               //                    cout<<"ijkcurr = "<<ijkcurr[0]<<" "<<ijkcurr[1]<<" "<<ijkcurr[2]<<endl;
               *(poctree->cartesian2octree(ijkcurr, levelMax)) = 1;
            }
         }
      }
   }

   //    octreeNarrowBand.printActivity();
   oKeyManager key_managerNB(*poctree);
   ExportGMSHAsciiNodes(*poctree, key_managerNB, "active_nodes_narrow_band_local");
}

void surface2LevelSetNarrowBand::activateAllCells() { poctree->activateFinestLevel(); }

void surface2LevelSetNarrowBand::computeLsOnNarrowBand(vector<double> *lsvVec)
{
   cout << "computeLsOnNarrowBand\n";

   // ANN - POINTS
   //    xgeom::xDistanceNearestPoint < xIter, xgeom::xAOMDToANN > distance_nearest_point(mesh_surf.begin(0), mesh_surf.end(0));
   // CGAL - TRIS
   xgeom::xDistanceNearestPoint<xIter, xgeom::xAOMDTriToCGALTri> distance_nearest_point(mesh_surf.begin(2), mesh_surf.end(2));
   // CGAL - EDGES
   //    xgeom::xDistanceNearestPoint < xIter, xgeom::xAOMDEdgeToCGALEdge > distance_nearest_point(mesh_surf.begin(1),
   //    mesh_surf.end(1));
   // CGAL - POINTS
   //    xgeom::xDistanceNearestPoint < xIter, xgeom::xAOMDVertexToCGALPoint > distance_nearest_point(mesh_surf.begin(0),
   //    mesh_surf.end(0));

   if (levelMax == -1) throw;

   // int dim=3;
   // int sizeMax = (int) ((pow(2.,levelMax)+1)* (pow(2.,levelMax)+1) * (pow(2.,levelMax)+1));

   int imax = (int)pow(2., levelMax) + 1;
   int jmax = (int)pow(2., levelMax) + 1;
   int kmax = (int)pow(2., levelMax) + 1;

   if (lsvVec) lsvVec->reserve(imax * jmax * kmax);

   // double hx=(BBmax(0)-BBmin(0))/(imax-1);
   // double hy=(BBmax(1)-BBmin(1))/(jmax-1);
   // double hz=(BBmax(2)-BBmin(2))/(kmax-1);

   double d;
   mEntity *nearestEntity;
   xtensor::xPoint closestPointOnsurf;

   bool verbose(0);

   // Compute distance for active nodes...
   oKeyManager key_managerNB(*poctree);
   oKeyManager::const_iterator it = key_managerNB.begin();
   oKeyManager::const_iterator ite = key_managerNB.end();
   if (lsvVec && static_cast<int>(lsvVec->size()) != imax * jmax * kmax) lsvVec->resize(imax * jmax * kmax);

   double P[3];
   for (; it != ite; ++it)
   {
      const oKey *key = *it;
      pmapping->ijk2xyz(key->getIJK(), levelMax, P);
      int nodeID = key->getId();

      if (verbose) cout << "gridpoint " << xPoint(P[0], P[1], P[2]) << endl;

      mVertex V(0, Trellis_Util::mPoint(P[0], P[1], P[2]), nullptr);
      distance_nearest_point.nearestEntityDistance(&V, nearestEntity, closestPointOnsurf, d);
      if (verbose) cout << "nearestpoint " << closestPointOnsurf << endl;

#if 1
      xtensor::xVector<> direction(closestPointOnsurf, V.point());
      std::pair<int, int> subEntity = getSubEntity(nearestEntity, closestPointOnsurf, verbose);
      if (verbose) cout << subEntity.first << " " << subEntity.second << endl;
      if (verbose) nearestEntity->print();

      xtensor::xVector<> normale;
      if (useWeightedPseudoNormals)
      {
         if (subEntity.first < 2)
            normale = weightedNormalsComp.getNormale(nearestEntity->get(subEntity.first, subEntity.second));
         else
            normale = weightedNormalsComp.getNormale(nearestEntity);

         if (normale.mag() < 1.e-12) verbose = 1;
         if (verbose) cout << "norm " << normale << " direction " << direction << endl;

         if (normale.mag() < 1.e-12)
         {
            // Cas pathologique de normale quasi nulle...
            // On prend la normale du segement ?
            if (verbose) cout << "vanishing normal... your surface may be non manifold\n";
            if (subEntity.first == 1) throw;  // Not coded for vanishing edge normals

            normale = weightedNormalsComp.getNormale(nearestEntity->get(subEntity.first + 1, subEntity.second));
            if (verbose) cout << "New norm " << normale << endl;
         }
      }
      else
      {
         normale = weightedNormalsComp.normaleToTri(nearestEntity);
      }
      double dir = direction * normale;
      if (verbose) cout << "dir " << dir << endl;
      if (dir < 0.)
         dir = -1.;
      else
         dir = 1.;

      (*lsvVec)[nodeID] = dir * d;
      verbose = 0;
#endif
   }
}

void surface2LevelSetNarrowBand::extendLset(vector<double> *lsvVec)
{
   cout << "extendLset\n";

   typedef unsigned char byte;
   std::vector<byte> activityVector(lsvVec->size(), 0);

   int imax = (int)pow(2., levelMax) + 1;
   int jmax = (int)pow(2., levelMax) + 1;
   int kmax = (int)pow(2., levelMax) + 1;
   int c[3] = {1, (topo.index_max_for_levels[0][levelMax] + 2),
               (topo.index_max_for_levels[0][levelMax] + 2) * (topo.index_max_for_levels[1][levelMax] + 2)};

   //   oKeyManager::container frontNodes;

   oKeyManager key_managerNB(*poctree);
   oKeyManager::const_iterator it = key_managerNB.begin();
   oKeyManager::const_iterator ite = key_managerNB.end();
   //  key_managerNB.printForDebug(cout);

   // set frozen nodes
   for (; it != ite; ++it)
   {
      activityVector[(*it)->getId()] = 1;
   }

   // Work on front
   it = key_managerNB.begin();
   //  oKeyManager::container frontNodes(it, ite);

   int neighbors[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
   // Travailler directement sur le keymanager ??
   // int iloop = 0;

   // ULTRA MOCHE : MODIFIER L'ALGO !!!
   // Eviter les hash (un rehash invalide l'iterateur du for(oKeyManager::const_iterator itf = k....)
   float LF = key_managerNB.max_load_factor();
   key_managerNB.max_load_factor(10 * LF);

   //  while(!frontNodes.empty()){
   while (key_managerNB.size())
   {
      //            cout<<"Marching "<<iloop++<<" "<<key_managerNB.size()<<endl;
      // Loop on front nodes
      for (oKeyManager::const_iterator itf = key_managerNB.begin(); itf != ite;)
      {
         //          cout<<key_managerNB.size()<<endl;
         const oKey *key = *itf;
         const int *ijk = key->getIJK();
         int iDf = key->getId();
         double valFront = (*lsvVec)[iDf];
         double valNeigh = 100.;
         if (valFront < 0.) valNeigh = -100.;

         // Loop Neighbor
         for (int iN = 0; iN < 6; ++iN)
         {
            int ijkN[3] = {ijk[0] + neighbors[iN][0], ijk[1] + neighbors[iN][1], ijk[2] + neighbors[iN][2]};
            if (ijkN[0] >= imax) ijkN[0] = imax - 1;
            if (ijkN[1] >= jmax) ijkN[1] = jmax - 1;
            if (ijkN[2] >= kmax) ijkN[2] = kmax - 1;

            if (ijkN[0] < 0) ijkN[0] = 0;
            if (ijkN[1] < 0) ijkN[1] = 0;
            if (ijkN[2] < 0) ijkN[2] = 0;

            int iDn = std::inner_product(ijkN, ijkN + 3, c, 0);
            //              cout<<iDn<<endl;

            // Si voisin inactif
            if (!activityVector[iDn])
            {
               (*lsvVec)[iDn] = valNeigh;
               key_managerNB.insert(ijkN);
               activityVector[iDn] = 1;
            }
         }

         ++itf;
         key_managerNB.erase(ijk);
      }
   }
}

void surface2LevelSetNarrowBand::saveLset(string lsname, vector<double> *lsvVec)
{
   cout << "saveLset\n";
   fstream file(lsname.c_str(), ios::in | ios::out | ios::trunc);

   int dim = 3;
   // int sizeMax = (int) ((pow(2.,levelMax)+1)* (pow(2.,levelMax)+1) * (pow(2.,levelMax)+1));

   int imax = (int)pow(2., levelMax) + 1;
   int jmax = (int)pow(2., levelMax) + 1;
   int kmax = (int)pow(2., levelMax) + 1;

   double hx = (BBmax(0) - BBmin(0)) / (imax - 1);
   double hy = (BBmax(1) - BBmin(1)) / (jmax - 1);
   double hz = (BBmax(2) - BBmin(2)) / (kmax - 1);

   file << dim << endl;       // dim
   file << imax << endl;      // sizeX
   file << jmax << endl;      // sizeY
   file << kmax << endl;      // sizeZ
   file << hx << endl;        // spaceX
   file << hy << endl;        // spaceY
   file << hz << endl;        // spaceZ
   file << BBmin(0) << endl;  // originX
   file << BBmin(1) << endl;  // originY
   file << BBmin(2) << endl;  // originZ

   // Beark...
   for (vector<double>::const_iterator i = lsvVec->begin(); i != lsvVec->end(); ++i)
   {
      file << *i << '\n';
   }
}

std::pair<int, int> surface2LevelSetNarrowBand::getSubEntity(mEntity *nearestEntity, const xtensor::xPoint &closestPointOnsurf,
                                                             bool verbose)
{
   xfem::xGeomElem geo(nearestEntity);
   geo.setUVWForXYZ(closestPointOnsurf);
   xtensor::xPoint uvw = geo.getUVW();

   double u(abs(uvw(0))), v(abs(uvw(1)));
   // double epsilon = 2. * std::numeric_limits<double>::epsilon();
   double epsilon = 1.e-13;
   if (verbose) cout << "epsilon=" << epsilon << endl;
   if (verbose) cout << "u=" << u << " v=" << v << endl;

   if (u >= epsilon && v >= epsilon && u + v <= 1. - epsilon)
   {
      return make_pair(2, 0);
   }
   else if (u < epsilon)
   {
      if (v < epsilon) return make_pair(0, 0);
      if (v > 1. - epsilon) return make_pair(0, 2);
      return make_pair(1, 2);
   }
   else if (v < epsilon)
   {
      if (u < epsilon) return make_pair(0, 0);
      if (u > 1. - epsilon) return make_pair(0, 1);
      return make_pair(1, 0);
   }
   else if (u + v > 1. - epsilon)
   {
      if (u < epsilon) return make_pair(0, 2);
      if (u > 1. - epsilon) return make_pair(0, 1);
      return make_pair(1, 1);
   }

   cout << "HaraKiri\n";
   throw;
}

#endif
#endif

}  // namespace xoctree
}  // namespace xinterface
