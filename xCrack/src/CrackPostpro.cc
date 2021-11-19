/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "CrackPostpro.h"

#include <iostream>
#include <map>
#include <vector>

#include "lCrack.h"
#include "mEntity.h"
#include "mTensor2.h"
#include "mTensor4.h"
#include "mVector.h"
#include "mVertex.h"
#include "xAlgorithm.h"
#include "xElement.h"
#include "xExportGmsh.h"
#include "xField.h"
#include "xForm.h"
#include "xGeomElem.h"
#include "xMaterial.h"
#include "xMaterialManager.h"
#include "xMesh.h"
#include "xSingleton.h"
#include "xTensors.h"
#include "xVectorField.h"
#include "xZone.h"

using namespace std;

using namespace xfem;
using namespace AOMD;

namespace xcrack
{
// pour modifier la fonction q
#define QMODIF 0
#define TCUBE (1. / 2.)

#define QRADIAL 0

// boucle sur les boîtes du front
// boucle sur les éléments de chaque boîte
// boucle sur les points d'intégration de chaque élt
// boucles sur les éléments finis touchant ces points

void lCrack::extractSIFs(const xField<>& disp_l, std::string filename, bool autop, double ratiop, double hp, double r1p,
                         double r2p)  // output file with the sifs
{
   const bool debug = false;
   //
   xEvalGradField<xtool::xIdentity<xtensor::xTensor2<>>> grad_disp(disp_l);
   char s2[10];
   sprintf(s2, "-%d", timestep);
   filename = filename + s2;
   filename = filename + ".txt";

   xexport::xExportGmshAscii pexport;

   // dimensions of the JDomain volume
   std::map<int, xtensor::xVector<>> JBoxDimension;  // for 3D
   // coordinates of the point on the crack front where J is computed
   std::map<int, xtensor::xPoint> JFrontPoint;
   // rotation matrix (X = Rx) from global x to local crack-front X coordinate system
   std::map<int, xtensor::xTensor2<>> JMap;

   xMesh* mesh_front1D = mesh_crack_front;
   xtensor::xTensor4<> DC;
   std::vector<double> MatInfo(2);
   //
   // Compute the length of each element on the front
   //
   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> length;
   // std::multimap<double,mEntity*> length_inv;
   // std::multimap<double,mEntity*>::value_type value_type;
   if (debug) cout << "mesh_front1D->begin(1) is " << mesh_front1D->size(1) << endl;
   for (xIter it = mesh_front1D->begin(1); it != mesh_front1D->end(1); ++it)
   {
      mEntity* e = *it;
      mVertex* v0 = (mVertex*)e->get(0, 0);
      xtensor::xPoint p0 = v0->point();
      mVertex* v1 = (mVertex*)e->get(0, 1);
      xtensor::xPoint p1 = v1->point();
      xtensor::xVector<> v(p0, p1);
      double d = std::sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2));
      length[e] = d;
      // std::pair<double, mEntity*> value(d,e);
      // length_inv.insert(value);
   }
   //
   // Compute the norm of the normal at each node on the 1D front
   //
   // const int was_created_by_tag = mesh->get_was_created_by_tag();
   std::unordered_map<mVertex*, double, EntityHashKey, EntityEqualKey> normal_norm;
   for (xIter it = mesh_front1D->begin(0); it != mesh_front1D->end(0); ++it)
   {
      mVertex* v = (mVertex*)*it;
      xtensor::xPoint p = v->point();
      mEntity* e2d = xMesh::get_const_was_created_by().at(*v);
      mEntity* e3d = xMesh::get_const_was_created_by().at(*e2d);
      if (e3d->getLevel() != 3) e3d = e3d->get(3, 0);
      xtensor::xVector<> e1, e2, e3;
      xtensor::xTensor2<> c1, c2, c3;
      xfem::xGeomElem geom(e3d);
      geom.setUVWForXYZ(p);
      getLocalCurv(e3d, geom.getUVW(), e1, e2, e3, c1, c2, c3);  // todo
      normal_norm[v] = std::sqrt(e1 * e1);
   }
   //
   // Compute the proper width for the domain integral
   //
   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> box_width;
   std::unordered_map<mEntity*, xtensor::xPoint, EntityHashKey, EntityEqualKey> box_location;
#if 0  // box around the front nodes
  std::vector<double> l(2);
  for (xIter it = mesh_front1D->begin(0); it != mesh_front1D->end(0); ++it ) {
    mVertex* v = (mVertex*) *it;
    if (v->size(1) == 2) {
      for (int j = 0; j < v->size(1); j++) {
	mEntity* e = v->get(1,j);
	l[j] = length[e];      
      }
      box_width[v] = std::min(l[0], l[1]);
      box_location[v] = v->point();
    }
  }
#endif
#if 1  // box around each edge of the 1D front
   std::vector<double> l(2);
   for (xIter it = mesh_front1D->begin(1); it != mesh_front1D->end(1); ++it)
   {
      mEntity* e = *it;
      box_width[e] = length.find(e)->second;
      mVertex* v0 = (mVertex*)e->get(0, 0);
      mVertex* v1 = (mVertex*)e->get(0, 1);
      box_location[e] = (v0->point() + v1->point()) * 0.5;
   }
#endif

   //  xVectorField speed_front(mesh_front1D);
   v1D.setSupport(mesh_front1D);

   xVectorField Ks(mesh_front1D);
   xLevelSet Js(mesh_front1D);
   xLevelSet angle(mesh_front1D);
   xVectorField position(mesh_front1D);

   // postid counts the crack front position
   int postid = 0;
   //  unsigned int dist = 0;
   bool out_of_the_box;
   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::iterator it = box_width.begin();
   std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey>::iterator itEnd = box_width.end();

   out_of_the_box = false;

   // printf("box_width size is %d\n", box_width.size());
   // loop over the boxes of the 1D front

   FILE* OUTPUT2 = fopen(filename.c_str(), "a");
   fprintf(OUTPUT2, "x y z r1 r2 r3 K1 K2 K3 J J(K1,K2,K3) K1*(J)\n");
   fclose(OUTPUT2);

   for (; it != itEnd; it++)
   {
      //  for ( ; it != itEnd;   it =  box_width.end()) {
      //       advance(it, min(box_width.size()/20, dist))) {
      //       advance(it, min(min(box_width.size()/9,1), dist))) {
      //    dist = 0;
      //    std::distance(it, itEnd, dist);

      mEntity* ebox = it->first;
      // initialize
      double Jhh = 0.0, vol = 0.0;
      xtensor::xVector<> Ih(0.);
      xtensor::xTensor2<> I(0.);

      // find the element the node is in to determine the size
      // of the box
      mEntity* e2d = xMesh::get_const_was_created_by().at(*ebox);
      mEntity* e3d = xMesh::get_const_was_created_by().at(*e2d);
      // find the average mesh size on the front
      double h = 0.0;
      // printf("Hellu3\n");
      if (e3d->getLevel() == 3)
      {
         xElement elem(e3d);
         h = pow(elem.getVolume(), 1. / 3.);
      }
      else
      {
         for (int j = 0; j < e3d->size(3); j++)
         {
            // printf("Hellu4\n");
            mEntity* e = e3d->get(3, j);
            // printf("Hellu5\n");
            xElement elem(e);
            h += pow(elem.getVolume(), 1. / 3.) / (double)e3d->size(3);
         }
      }
      if (debug) cout << "mesh size near the front h " << h << endl;

      // h = 2.5e-02;
      // h *= 1.5;
      // double r3 = 2 * h;
      // tests

      double r3, r1, r2;
      if (autop)
      {
         r3 = h;
         r1 = ratiop * h;
         r2 = ratiop * h;
      }
      else
      {
         r3 = hp;
         r1 = r1p;
         r2 = r2p;
      }

      cout << "the box size  is " << r1 << " " << r2 << " " << r3 << endl;

      // creating the elements we will integrate on
      std::map<int, xtensor::xPoint> JCoord;
      std::map<int, std::vector<int>> JElements;
      xtensor::xPoint p = box_location.find(ebox)->second;
      if (e3d->getLevel() != 3) e3d = e3d->get(3, 0);
      xtensor::xPoint local;
      xtensor::xVector<> eo1, eo2, eo3;
      xElement elem(e3d);
      elem.xyz2uvw(p);
      getLocalSmoothOrthoAxis(e3d, elem.getUvw(), local, eo1, eo2, eo3);  // todo

      xtensor::xTensor2<> axis;
      for (int j = 0; j < 3; j++)
      {
         axis(0, j) = eo1(j);
         axis(1, j) = eo2(j);
         axis(2, j) = eo3(j);
      }
      // double pointJfront[3];
      xtensor::xPoint pointJfront(p);  //(0); pointJfront[1] = p(1); pointJfront[2] = p(2);
      xtensor::xTensor2<> transaxis;
      // matTrans(axis,transaxis);
      transaxis = axis;
      Transpo(transaxis);
      JGrid3D(pointJfront, r1, r2, r3, transaxis, JCoord, JElements);

      xMesh* postpro_mesh = new xMesh;

      CreateMeshWithPostProElements(postpro_mesh, JCoord, JElements);

      // data->SetJBox(postid, r1, r2, r3);
      JBoxDimension.insert(make_pair(postid, xtensor::xVector<>(r1, r2, r3)));

      // data->SetJFrontPoint(postid, pointJfront);
      JFrontPoint.insert(make_pair(postid, pointJfront));

      // data->SetJMap(postid, axis);
      JMap.insert(make_pair(postid, axis));
      if (debug) cout << "size of postpro_mesh is " << postpro_mesh->size(postpro_mesh->dim()) << endl;
      // loop over the elements for the domain integral
      for (xIter itj = postpro_mesh->begin(postpro_mesh->dim()); itj != postpro_mesh->end(postpro_mesh->dim()); ++itj)
      {
         mEntity* ej = *itj;
         xfem::xGeomElem JEM(ej);
         // it faut retrouver 216 points JEM.SetGaussPointNumber(216);
         JEM.SetIntegrationPointNumberForDegree(10);  // 10 is the right degree for 216 gauss points
         if (debug)
         {
            printf("No of quadrature points is %d \n", JEM.GetNbIntegrationPoints());
            printf("Weight and DetJac are: %22.15e %2.15e \n", JEM.GetWeight(), JEM.GetDetJac());
         }

         for (int i = 0; i < static_cast<int>(JEM.GetNbIntegrationPoints()); ++i)
         {
            JEM.setUVW(i);  // Upos, Vpos, Wpos and JacMatrix are set
            xtensor::xPoint xyz = JEM.getXYZ();
            if (debug) printf("XPos YPos ZPos %12.5e %12.5e %12.5e\n", xyz(0), xyz(1), xyz(2));
            auto list_pe_uvw = mesh->locateElement(xyz);  // todo

            if (list_pe_uvw.empty()) out_of_the_box = true;

            double eset_size = (double)std::distance(list_pe_uvw.begin(), list_pe_uvw.end());
            for (auto pe_uvw : list_pe_uvw)
            {
               const AOMD::mEntity* pe = pe_uvw.first;
               xfem::xGeomElem geom_appro(const_cast<AOMD::mEntity*>(pe));
               xMaterial* mat = xMaterialManagerSingleton::instance().getMaterial(&geom_appro);
               const xTensors* properties = mat->getProperties();
               MatInfo[0] = properties->scalar("YOUNG_MODULUS");
               MatInfo[1] = properties->scalar("POISSON_RATIO");
               mat->sensitivityTo("strain", DC);

               xtensor::xPoint xyzloc;
               xtensor::xTensor2<> matgradqloc;
               GetGradqLoc(postid, xyz, matgradqloc, xyzloc, JBoxDimension, JFrontPoint, JMap);
               xtensor::xVector<> gradqloc, dalpha;
               for (int j = 0; j < 3; j++) gradqloc(j) = matgradqloc(0, j);
               xtensor::xTensor2<> axis = JMap.find(postid)->second;
               dalpha = gradqloc * axis;
               double alpha = GetqLoc(postid, xyz, JBoxDimension, JFrontPoint, JMap);
               if (debug) cout << "alpha is " << alpha << endl;
               geom_appro.setUVW(pe_uvw.second);

               double Jhh_l;
               xtensor::xVector<> Ih_l;
               xtensor::xTensor2<> I_l;
               double vol_l;
               if (debug) cout << " entering GetJint3d" << endl;
               GetJint3DLevelset(MatInfo, grad_disp, alpha, dalpha, DC, JEM.GetWeight(), JEM.GetDetJac(), &geom_appro, this,
                                 Jhh_l, Ih_l, I_l, vol_l);
               if (debug) cout << " leaving GetJint3d" << endl;
               Jhh += Jhh_l / eset_size;
               I_l /= eset_size;
               I += I_l;
               Ih += (Ih_l / eset_size);
               vol += vol_l;
               if (debug)
               {
                  cout << "Jhh " << Jhh << endl;
                  cout << "Ih  " << Ih << endl;
                  cout << "I   " << endl;
                  cout << I << endl;
                  cout << "vol " << vol << endl;
               }
            }
         }
      }
      delete postpro_mesh;
      if (debug) cout << "volume of the box is " << vol << endl;
      //
      xtensor::xVector<> boxdim = JBoxDimension.find(postid)->second;
#if QMODIF == 3
      double denominator = 1.0 * boxdim(2);  // 1.0*L_3
#elif QMODIF == 2
      double denominator = (1.0 + TCUBE) / 2 * boxdim(2);  // *L_3
#elif QMODIF == 1
      double denominator = 1.0 * boxdim(2);  // 1.0*L_3
#else
      double denominator = 0.5 * boxdim(2);  // 0.5*L_3
#endif
      Jhh /= denominator;
      Ih /= denominator;
      I /= denominator;

      double Estar = MatInfo[0] / (1. - MatInfo[1] * MatInfo[1]);
      double mu = MatInfo[0] / (2. * (1. + MatInfo[1]));

      //    FILE* OUTPUT = fopen(filename.c_str(),"a");

      printf("***********************************************************************\n");
      if (out_of_the_box)
      {
         printf("point hors de la boite : valeurs surement fausses %d \n", postid);
      }

      printf("Crack is 1, Front is 1 and PostId is %d \n", postid);
      printf("Point %d on crack front: %15.7e %15.7e %15.7e\n", postid, pointJfront(0), pointJfront(1), pointJfront(2));
      printf("MatInfo is: %15.7e %15.7e\n", MatInfo[0], MatInfo[1]);

      printf("J-INTEGRAL RESULT: Jdom = %22.14e\n", Jhh);
      double K1_modeI = sqrt(Estar * Jhh);
      printf("K1 from J assuming pure mode I = %22.15e\n", K1_modeI);

      double K1aux = sqrt(Estar * I(0, 0));
      double K2aux = sqrt(Estar * I(1, 1));
      double K3aux = sqrt(2. * mu * I(2, 2));

      printf("Interaction integral results\n");
      double K1 = Estar * Ih(0) / (2.0);
      double K2 = Estar * Ih(1) / (2.0);
      double K3 = 2. * mu * Ih(2) / (2.0);
      printf("K1  = %22.15e\n", K1);
      printf("K2  = %22.15e\n", K2);
      printf("K3  = %22.15e\n", K3);
      double J_fromKs = (K1 * K1 + K2 * K2) / Estar + K3 * K3 / (2.0 * mu);
      printf("(K1*K1 + K2*K2)/Estar + K3*K3/(2.0*mu)   = %22.15e\n", J_fromKs);
      printf("Without numerical errors, the value above should equal Jdom\n");
      printf("K1aux  = %22.15e\n", K1aux);
      printf("K2aux  = %22.15e\n", K2aux);
      printf("K3aux  = %22.15e\n", K3aux);
      printf("I11  = %22.15e\n", I(0, 0));
      printf("I22  = %22.15e\n", I(1, 1));
      printf("I33  = %22.15e\n", I(2, 2));
      //      printf( "The table below is\n");
      //      printf( "I11 I12 I13\n");
      //      printf( "I12 I22 I23\n");
      //      printf( "I13 I23 I33\n");
      //      printf( "It should be close to identity\n");
      //      printf( "%22.15e %22.15e %22.15e\n", I(0,0) * Estar, I(0,1), I(0,2));
      //      printf( "%22.15e %22.15e %22.15e\n", I(1,0), I(1,1) * Estar, I(1,2));
      //      printf( "%22.15e %22.15e %22.15e\n", I(2,0), I(2,1), I(2,2)*2.*mu);

      //    fclose(OUTPUT);

      FILE* OUTPUT2 = fopen(filename.c_str(), "a");
      cout << "out_of the box" << out_of_the_box << endl;
      if (out_of_the_box == false)
      {
         fprintf(OUTPUT2, "%e %e %e %e %e %e %e %e %e %e %e %e\n", pointJfront(0), pointJfront(1), pointJfront(2), r1, r2, r3, K1,
                 K2, K3, Jhh, J_fromKs, K1_modeI);
      }
      out_of_the_box = false;
      fclose(OUTPUT2);

      postid++;

      // printf("before computing speed_front\n");
      // printf("eo1 is %12.5e %12.5e %12.5e\n", eo1(0), eo1(1), eo1(2));
      // printf("Jhh is %12.5e\n", Jhh);

      double theta = AngleWithMaxHoopStress(K1, K2);

      v1D(ebox) = (eo1 * cos(theta) + eo2 * sin(theta)) * Jhh;

      // speed_front(ebox) = (eo1 * cos(theta) + eo2 * sin(theta)) * J_fromKs;
      Ks(ebox) = xtensor::xVector<>(K1, K2, K3);
      position(ebox) = xtensor::xVector<>(p(0), p(1), p(2));
      const double PI = 4. * atan(1.);
      // printf("speed vector is %12.5e %12.5e %12.5e\n",
      //       speed_front(ebox)(0), speed_front(ebox)(1), speed_front(ebox)(2));
      Js(ebox) = Jhh;
      angle(ebox) = 180. * atan2(p(1), p(0)) / PI;
      printf("box size  %12.5e %12.5e %12.5e K1: %lf K2: %lf K3: %lf\n", r1, r2, r3, K1, K2, K3);
   }

   mesh_front1D->getMesh().modifyState(0, 1, true);
   for (xIter it = mesh_front1D->begin(0); it != mesh_front1D->end(0); ++it)
   {
      // printf("looping over the modes for the export\n");
      mVertex* v = (mVertex*)*it;
      v1D(v) = xtensor::xVector<>(0.0, 0.0, 0.0);
      for (int j = 0; j < v->size(1); j++)
      {
         mEntity* e = v->get(1, j);
         // we make the average at the node of the two neighbors speeds
         v1D(v) += (v1D(e) / (double)v->size(1));
      }
      // printf("speed vector is %12.5e %12.5e %12.5e\n",
      //       speed_front(v)(0), speed_front(v)(1), speed_front(v)(2));
   }
   // printf("before export\n");
   // toputback
   {
      string s = "speedfront";
      char s2[10];
      sprintf(s2, "-%d", timestep);
      s = s + s2;

      v1D.exportGmsh(s);
      // Export(v1D, pexport, s);
   }

   //    const double PI = 4.*atan(1.);
   //    for (mIteratorPtr it(mesh_front1D->getLevelIterator(1)); !it->done(); it->next() ) {
   //      mEntity* e = it->curr();
   //      xtensor::xPoint p = box_location.find(e)->second;
   //      angle(e) = 180. * atan2(p(1), p(0))/PI;
   //    }
   // 1 indicates we export on the edges
   /*
     std::ofstream fout("Ks.m");
     angle.exportMatlab(fout, "angle", 1);
     position.exportMatlab(fout, "position", 1);
     Ks.exportMatlab(fout, "K", 1);
     Js.exportMatlab(fout, "J", 1);
     fout.close();
   */
   // todo  crack.setFrontVelocity(speed_front);

   return;
}

// namespace xcrack
//{

xtensor::xTensor2<>& Transpo(xtensor::xTensor2<>& in)
{
   double temp = in(0, 1);
   in(0, 1) = in(1, 0);
   in(1, 0) = temp;
   temp = in(0, 2);
   in(0, 2) = in(2, 0);
   in(2, 0) = temp;
   temp = in(1, 2);
   in(1, 2) = in(2, 1);
   in(2, 1) = temp;
   return in;
}

void GetMech3DAuxFields(const xtensor::xPoint& XLoc, const vector<double>& MatInfo, xtensor::xTensor2<>& AuxStress,
                        xtensor::xTensor2<>& AuxEps, xtensor::xVector<>& AuxDisp, xtensor::xVector<>& AuxGrad0Disp,
                        xtensor::xVector<>& AuxGrad1Disp, xtensor::xTensor2<>& AuxGrad0Stress,
                        xtensor::xTensor2<>& AuxGrad1Stress, xtensor::xVector<>& AuxGrad00Disp, xtensor::xVector<>& AuxGrad11Disp,
                        xtensor::xVector<>& AuxGrad01Disp, int mode)
{
   //  - extension au 3D - prise en compte du mode III          - 12/2000
   //  - expressions de u1,du1dr,du1dt,u2,du2dr,du2dt
   //    simplifiees pour les modes I et II (Bui 1978)

   const bool debug = false;
   if (debug) cout << "In Auxiliary fields : Xloc " << XLoc << endl;
   const double PI = 4. * atan(1.);
   double r = sqrt((XLoc(0) * XLoc(0)) + (XLoc(1) * XLoc(1)));
   double theta = atan2(XLoc(1), XLoc(0));

   double young, nu, mu;
   int i, j;

   double CT, ST, CT2, ST2, C3T2, S3T2;
   double drdx, drdy, dtdx, dtdy;
   double d2rdxx, d2rdxy, d2rdyy, d2tdxx, d2tdxy, d2tdyy;
   double du1dr, du1dt, du2dr, du2dt, du3dr, du3dt;
   double ds11dr, ds11dt, ds12dr, ds12dt, ds22dr, ds22dt;
   double ds13dr, ds13dt, ds23dr, ds23dt;
   double d2u3drdr, d2u3drdt, d2u3dtdr, d2u3dtdt;
   double d2u1drdr, d2u1drdt, d2u1dtdr, d2u1dtdt, d2u2drdr, d2u2drdt, d2u2dtdr, d2u2dtdt;
   double K1, K2, K3;
   double FACStress1, FACDisp1, FACStress2, FACDisp2, FACStress3, FACDisp3;

   // for checking the constitutive relations
   double kappa;

   K1 = 1.0;
   K2 = 1.0;
   K3 = 1.0;

   // zero all auxiliary fields
   for (i = 0; i < 3; i++)
   {
      AuxDisp(i) = 0.;
      AuxGrad0Disp(i) = 0.;
      AuxGrad1Disp(i) = 0.;
      AuxGrad00Disp(i) = 0.;
      AuxGrad11Disp(i) = 0.;
      AuxGrad01Disp(i) = 0.;
      for (j = 0; j < 3; j++)
      {
         AuxStress(i, j) = 0.;
         AuxEps(i, j) = 0.;
         AuxGrad0Stress(i, j) = 0.;
         AuxGrad1Stress(i, j) = 0.;
      }
   }

   young = MatInfo[0];
   nu = MatInfo[1];
   if (debug) cout << "young " << young << " nu " << nu << endl;
   mu = young / (2. * (1. + nu));

   // DC = young/((1.+nu)*(1. - 2.*nu));

   CT = cos(theta);
   ST = sin(theta);
   CT2 = cos(theta / 2.0);
   ST2 = sin(theta / 2.0);
   C3T2 = cos(3 * theta / 2.0);
   S3T2 = sin(3 * theta / 2.0);

   //: deriveres premieres utiles
   drdx = cos(theta);
   drdy = sin(theta);
   dtdx = -sin(theta) / r;
   dtdy = cos(theta) / r;

   //: derivees secondes utiles (ca devient chaud...)
   d2rdxx = sin(theta) * sin(theta) / r;
   d2rdxy = -sin(theta) * cos(theta) / r;
   d2rdyy = cos(theta) * cos(theta) / r;
   d2tdxx = 2. * sin(theta) * cos(theta) / (r * r);
   d2tdxy = -1. / (r * r) + 2. * sin(theta) * sin(theta) / (r * r);
   d2tdyy = -2. * sin(theta) * cos(theta) / (r * r);

   FACStress1 = sqrt(1. / (2. * PI));
   FACDisp1 = sqrt(1. / (2. * PI)) / (2 * mu);
   FACStress2 = sqrt(1. / (2. * PI));
   FACDisp2 = sqrt(1. / (2. * PI)) / (2 * mu);
   FACStress3 = sqrt(1. / (2. * PI));
   FACDisp3 = 2. * sqrt(1. / (2. * PI)) / (mu);

   kappa = 3.0 - 4. * nu;

   switch (mode)
   {
      case 1:  // Auxiliary Fields - Mode I
         AuxStress(0, 0) = K1 * FACStress1 / sqrt(r) * CT2 * (1.0 - ST2 * S3T2);
         AuxStress(1, 1) = K1 * FACStress1 / sqrt(r) * CT2 * (1.0 + ST2 * S3T2);
         AuxStress(0, 1) = K1 * FACStress1 / sqrt(r) * ST2 * CT2 * C3T2;
         AuxStress(1, 0) = AuxStress(0, 1);

         // 16 avril 2001
         //  ... gradient des contraintes
         ds11dr = -K1 * FACStress1 / (2. * pow(r, 3. / 2.)) * CT2 * (1.0 - ST2 * S3T2);
         ds11dt = K1 * FACStress1 / (-2. * sqrt(r)) * (ST2 * (1.0 - ST2 * S3T2) + CT2 * (CT2 * S3T2 + 3. * ST2 * C3T2));
         ds22dr = -K1 * FACStress1 / (2. * pow(r, 3. / 2.)) * CT2 * (1.0 + ST2 * S3T2);
         ds22dt = K1 * FACStress1 / (-2. * sqrt(r)) * (ST2 * (1.0 + ST2 * S3T2) - CT2 * (CT2 * S3T2 + 3. * ST2 * C3T2));
         ds12dr = -K1 * FACStress1 / (2. * pow(r, 3. / 2.)) * CT2 * (ST2 * C3T2);
         ds12dt = K1 * FACStress1 / (-2. * sqrt(r)) * (ST2 * (ST2 * C3T2) + CT2 * (-CT2 * C3T2 + 3. * ST2 * S3T2));
         //  ... dans le repere local
         AuxGrad0Stress(0, 0) = ds11dr * drdx + ds11dt * dtdx;
         AuxGrad1Stress(0, 0) = ds11dr * drdy + ds11dt * dtdy;
         AuxGrad0Stress(0, 1) = ds12dr * drdx + ds12dt * dtdx;
         AuxGrad1Stress(0, 1) = ds12dr * drdy + ds12dt * dtdy;
         AuxGrad0Stress(1, 0) = AuxGrad0Stress(0, 1);
         AuxGrad1Stress(1, 0) = AuxGrad1Stress(0, 1);
         AuxGrad0Stress(1, 1) = ds22dr * drdx + ds22dt * dtdx;
         AuxGrad1Stress(1, 1) = ds22dr * drdy + ds22dt * dtdy;

         //
         AuxDisp(0) = K1 * FACDisp1 * sqrt(r) * CT2 * (kappa - CT);
         AuxDisp(1) = K1 * FACDisp1 * sqrt(r) * ST2 * (kappa - CT);

         du1dr = K1 * FACDisp1 * 0.5 / sqrt(r) * CT2 * (kappa - CT);
         du1dt = K1 * FACDisp1 * sqrt(r) * (-0.5 * ST2 * (kappa - CT) + CT2 * ST);
         du2dr = K1 * FACDisp1 * 0.5 / sqrt(r) * ST2 * (kappa - CT);
         du2dt = K1 * FACDisp1 * sqrt(r) * (0.5 * CT2 * (kappa - CT) + ST2 * ST);

         AuxGrad0Disp(0) = du1dr * drdx + du1dt * dtdx;
         AuxGrad1Disp(0) = du1dr * drdy + du1dt * dtdy;
         AuxGrad0Disp(1) = du2dr * drdx + du2dt * dtdx;
         AuxGrad1Disp(1) = du2dr * drdy + du2dt * dtdy;

         AuxEps(0, 0) = AuxGrad0Disp(0);
         AuxEps(1, 0) = 0.5 * (AuxGrad1Disp(0) + AuxGrad0Disp(1));
         AuxEps(0, 1) = AuxEps(1, 0);
         AuxEps(1, 1) = AuxGrad1Disp(1);

         // 16 avril 2001
         //  ... gradient du gradient des deplacements
         d2u1drdr = -K1 * FACDisp1 / 4. / pow(r, 3. / 2.) * CT2 * (kappa - CT);
         d2u1drdt = K1 * FACDisp1 / 2. / sqrt(r) * ((-0.5 * ST2) * (kappa - CT) + CT2 * ST);
         d2u1dtdr = d2u1drdt;
         d2u1dtdt = K1 * FACDisp1 * sqrt(r) * ((-0.25 * CT2) * (kappa - CT) + C3T2);
         d2u2drdr = -K1 * FACDisp1 / 4. / pow(r, 3. / 2.) * ST2 * (kappa - CT);
         d2u2drdt = K1 * FACDisp1 / 2. / sqrt(r) * ((0.5 * CT2) * (kappa - CT) + ST2 * ST);
         d2u2dtdr = d2u2drdt;
         d2u2dtdt = K1 * FACDisp1 * sqrt(r) * ((-0.25 * ST2) * (kappa - CT) + S3T2);
         //  ... dans le repere local
         AuxGrad00Disp(0) = (d2u1drdr * drdx + d2u1drdt * dtdx) * drdx + (d2u1dtdr * drdx + d2u1dtdt * dtdx) * dtdx +
                            (du1dr * d2rdxx + du1dt * d2tdxx);
         AuxGrad01Disp(0) = (d2u1drdr * drdx + d2u1drdt * dtdx) * drdy + (d2u1dtdr * drdx + d2u1dtdt * dtdx) * dtdy +
                            (du1dr * d2rdxy + du1dt * d2tdxy);
         AuxGrad11Disp(0) = (d2u1drdr * drdy + d2u1drdt * dtdy) * drdy + (d2u1dtdr * drdy + d2u1dtdt * dtdy) * dtdy +
                            (du1dr * d2rdyy + du1dt * d2tdyy);
         AuxGrad00Disp(1) = (d2u2drdr * drdx + d2u2drdt * dtdx) * drdx + (d2u2dtdr * drdx + d2u2dtdt * dtdx) * dtdx +
                            (du2dr * d2rdxx + du2dt * d2tdxx);
         AuxGrad01Disp(1) = (d2u2drdr * drdx + d2u2drdt * dtdx) * drdy + (d2u2dtdr * drdx + d2u2dtdt * dtdx) * dtdy +
                            (du2dr * d2rdxy + du2dt * d2tdxy);
         AuxGrad11Disp(1) = (d2u2drdr * drdy + d2u2drdt * dtdy) * drdy + (d2u2dtdr * drdy + d2u2dtdt * dtdy) * dtdy +
                            (du2dr * d2rdyy + du2dt * d2tdyy);

         break;
      case 2:  // Auxiliary Fields - Mode II
         AuxStress(0, 0) = -K2 * FACStress2 / sqrt(r) * ST2 * (2.0 + CT2 * C3T2);
         AuxStress(1, 1) = K2 * FACStress2 / sqrt(r) * ST2 * CT2 * C3T2;
         AuxStress(0, 1) = K2 * FACStress2 / sqrt(r) * CT2 * (1.0 - ST2 * S3T2);
         AuxStress(1, 0) = AuxStress(0, 1);

         // 16 avril 2001
         //  ... gradient des contraintes
         ds11dr = K2 * FACStress2 / (2. * pow(r, 3. / 2.)) * ST2 * (2.0 + CT2 * C3T2);
         ds11dt = -K2 * FACStress2 / (2. * sqrt(r)) * (CT2 * (2.0 + CT2 * C3T2) - ST2 * (ST2 * C3T2 + 3. * CT2 * S3T2));
         ds22dr = -K2 * FACStress2 / (2. * pow(r, 3. / 2.)) * ST2 * (CT2 * C3T2);
         ds22dt = K2 * FACStress2 / (2. * sqrt(r)) * (CT2 * (CT2 * C3T2) - ST2 * (ST2 * C3T2 + 3. * CT2 * S3T2));
         ds12dr = -K2 * FACStress2 / (2. * pow(r, 3. / 2.)) * CT2 * (1.0 - ST2 * S3T2);
         ds12dt = K2 * FACStress2 / (2. * sqrt(r)) * (-ST2 * (1.0 - ST2 * S3T2) - CT2 * (CT2 * S3T2 + 3. * ST2 * C3T2));
         //  ... dans le repere local
         AuxGrad0Stress(0, 0) = ds11dr * drdx + ds11dt * dtdx;
         AuxGrad1Stress(0, 0) = ds11dr * drdy + ds11dt * dtdy;
         AuxGrad0Stress(0, 1) = ds12dr * drdx + ds12dt * dtdx;
         AuxGrad1Stress(0, 1) = ds12dr * drdy + ds12dt * dtdy;
         AuxGrad0Stress(1, 0) = AuxGrad0Stress(0, 1);
         AuxGrad1Stress(1, 0) = AuxGrad1Stress(0, 1);
         AuxGrad0Stress(1, 1) = ds22dr * drdx + ds22dt * dtdx;
         AuxGrad1Stress(1, 1) = ds22dr * drdy + ds22dt * dtdy;

         //
         AuxDisp(0) = K2 * FACDisp2 * sqrt(r) * ST2 * (kappa + 2 + CT);
         AuxDisp(1) = -K2 * FACDisp2 * sqrt(r) * CT2 * (kappa - 2 + CT);

         du1dr = K2 * FACDisp2 * 0.5 / sqrt(r) * ST2 * (kappa + 2 + CT);
         du1dt = K2 * FACDisp2 * sqrt(r) * (0.5 * CT2 * (kappa + 2 + CT) - ST2 * ST);
         du2dr = -K2 * FACDisp2 * 0.5 * (1. / sqrt(r)) * CT2 * (kappa - 2 + CT);
         du2dt = -K2 * FACDisp2 * sqrt(r) * (-0.5 * ST2 * (kappa - 2 + CT) - CT2 * ST);

         AuxGrad0Disp(0) = du1dr * drdx + du1dt * dtdx;
         AuxGrad1Disp(0) = du1dr * drdy + du1dt * dtdy;
         AuxGrad0Disp(1) = du2dr * drdx + du2dt * dtdx;
         AuxGrad1Disp(1) = du2dr * drdy + du2dt * dtdy;

         AuxEps(0, 0) = AuxGrad0Disp(0);
         AuxEps(1, 0) = 0.5 * (AuxGrad1Disp(0) + AuxGrad0Disp(1));
         AuxEps(0, 1) = AuxEps(1, 0);
         AuxEps(1, 1) = AuxGrad1Disp(1);

         // 16 avril 2001
         //  ... gradient du gradient des deplacements
         d2u1drdr = -K2 * FACDisp2 / 4. / pow(r, 3. / 2.) * ST2 * (kappa + 2.0 + CT);
         d2u1drdt = K2 * FACDisp2 / 2. / sqrt(r) * ((0.5 * CT2) * (kappa + 2.0 + CT) - ST2 * ST);
         d2u1dtdr = d2u1drdt;
         d2u1dtdt = K2 * FACDisp2 * sqrt(r) * ((-0.25 * ST2) * (kappa + 2.0 + CT) - S3T2);
         d2u2drdr = K2 * FACDisp2 / 4. / pow(r, 3. / 2.) * CT2 * (kappa - 2.0 + CT);
         d2u2drdt = -K2 * FACDisp2 / 2. / sqrt(r) * ((-0.5 * ST2) * (kappa - 2.0 + CT) - CT2 * ST);
         d2u2dtdr = d2u2drdt;
         d2u2dtdt = -K2 * FACDisp2 * sqrt(r) * ((-0.25 * CT2) * (kappa - 2.0 + CT) - C3T2);
         //  ... dans le repere local
         AuxGrad00Disp(0) = (d2u1drdr * drdx + d2u1drdt * dtdx) * drdx + (d2u1dtdr * drdx + d2u1dtdt * dtdx) * dtdx +
                            (du1dr * d2rdxx + du1dt * d2tdxx);
         AuxGrad01Disp(0) = (d2u1drdr * drdx + d2u1drdt * dtdx) * drdy + (d2u1dtdr * drdx + d2u1dtdt * dtdx) * dtdy +
                            (du1dr * d2rdxy + du1dt * d2tdxy);
         AuxGrad11Disp(0) = (d2u1drdr * drdy + d2u1drdt * dtdy) * drdy + (d2u1dtdr * drdy + d2u1dtdt * dtdy) * dtdy +
                            (du1dr * d2rdyy + du1dt * d2tdyy);
         AuxGrad00Disp(1) = (d2u2drdr * drdx + d2u2drdt * dtdx) * drdx + (d2u2dtdr * drdx + d2u2dtdt * dtdx) * dtdx +
                            (du2dr * d2rdxx + du2dt * d2tdxx);
         AuxGrad01Disp(1) = (d2u2drdr * drdx + d2u2drdt * dtdx) * drdy + (d2u2dtdr * drdx + d2u2dtdt * dtdx) * dtdy +
                            (du2dr * d2rdxy + du2dt * d2tdxy);
         AuxGrad11Disp(1) = (d2u2drdr * drdy + d2u2drdt * dtdy) * drdy + (d2u2dtdr * drdy + d2u2dtdt * dtdy) * dtdy +
                            (du2dr * d2rdyy + du2dt * d2tdyy);

         break;
      case 3:  // Auxiliary Fields - Mode III
         AuxStress(0, 2) = -K3 * FACStress3 / sqrt(r) * ST2;
         AuxStress(1, 2) = K3 * FACStress3 / sqrt(r) * CT2;
         AuxStress(2, 0) = AuxStress(0, 2);
         AuxStress(2, 1) = AuxStress(1, 2);

         // 16 avril 2001
         //  ... gradient des contraintes
         ds13dr = K3 * FACStress3 / 2. / pow(r, 3. / 2.) * ST2;
         ds13dt = -K3 * FACStress3 / sqrt(r) * (0.5) * CT2;
         ds23dr = -K3 * FACStress3 / 2. / pow(r, 3. / 2.) * CT2;
         ds23dt = -K3 * FACStress3 / sqrt(r) * (0.5) * ST2;
         //  ... dans le repere local
         AuxGrad0Stress(0, 2) = ds13dr * drdx + ds13dt * dtdx;
         AuxGrad1Stress(0, 2) = ds13dr * drdy + ds13dt * dtdy;
         AuxGrad0Stress(1, 2) = ds23dr * drdx + ds23dt * dtdx;
         AuxGrad1Stress(1, 2) = ds23dr * drdy + ds23dt * dtdy;
         AuxGrad0Stress(2, 1) = AuxGrad0Stress(1, 2);
         AuxGrad1Stress(2, 1) = AuxGrad1Stress(1, 2);
         AuxGrad0Stress(2, 0) = AuxGrad0Stress(0, 2);
         AuxGrad1Stress(2, 0) = AuxGrad1Stress(0, 2);

         //
         AuxDisp(2) = K3 * FACDisp3 * sqrt(r) * ST2;

         du3dr = K3 * FACDisp3 * 0.5 / sqrt(r) * ST2;
         du3dt = K3 * FACDisp3 * sqrt(r) * (0.5 * CT2);

         AuxGrad0Disp(2) = du3dr * drdx + du3dt * dtdx;
         AuxGrad1Disp(2) = du3dr * drdy + du3dt * dtdy;

         AuxEps(0, 2) = 0.5 * (0.0 + AuxGrad0Disp(2));
         AuxEps(1, 2) = 0.5 * (0.0 + AuxGrad1Disp(2));
         AuxEps(2, 0) = AuxEps(0, 2);
         AuxEps(2, 1) = AuxEps(1, 2);

         // 16 avril 2001
         //  ... gradient du gradient des deplacements
         d2u3drdr = -K3 * FACDisp3 / (4. * pow(r, 3. / 2.)) * ST2;
         d2u3drdt = K3 * FACDisp3 / (4. * sqrt(r)) * CT2;
         d2u3dtdr = d2u3drdt;
         d2u3dtdt = -K3 * FACDisp3 / 4. * sqrt(r) * ST2;
         //  ... dans le repere local
         AuxGrad00Disp(2) = (d2u3drdr * drdx + d2u3drdt * dtdx) * drdx + (d2u3drdt * drdx + d2u3dtdt * dtdx) * dtdx +
                            (du3dr * d2rdxx + du3dt * d2tdxx);
         AuxGrad01Disp(2) = (d2u3drdr * drdx + d2u3dtdr * dtdx) * drdy + (d2u3dtdr * drdx + d2u3dtdt * dtdx) * dtdy +
                            (du3dr * d2rdxy + du3dt * d2tdxy);
         AuxGrad11Disp(2) = (d2u3drdr * drdy + d2u3drdt * dtdy) * drdy + (d2u3drdt * drdy + d2u3dtdt * dtdy) * dtdy +
                            (du3dr * d2rdyy + du3dt * d2tdyy);

         break;
      default:
         assert(1 == 0);
         break;
   }

   // sigma33
   AuxStress(2, 2) = nu * (AuxStress(0, 0) + AuxStress(1, 1));
   AuxGrad0Stress(2, 2) = nu * (AuxGrad0Stress(0, 0) + AuxGrad0Stress(1, 1));
   AuxGrad1Stress(2, 2) = nu * (AuxGrad1Stress(0, 0) + AuxGrad1Stress(1, 1));

   return;
}

void GetJint3DLevelset(const std::vector<double>& MatInfo, xEval<xtensor::xTensor2<>>& grad_disp, double alpha,
                       const xtensor::xVector<>& dalpha2, const xtensor::xTensor4<>& Phys, double Weight, double DetJac,
                       xfem::xGeomElem* geo_elem, const lCrack* LS_Crack, double& Jhh, xtensor::xVector<>& Ih,
                       xtensor::xTensor2<>& I, double& vol)
{
   const bool debug = false;
   if (debug) cout << "DetJac in  GetJint3DLevelset " << DetJac << endl;
   // geometrical information at the point of integration
   //
   LS_Crack->setLocalInfos(geo_elem->getEntity(), geo_elem->getUVW());
   /// auxiliaryFields->SetLocation(LS_Crack->local);
   xtensor::xVector<>& no1 = LS_Crack->no1;
   //  xtensor::xVector<>&  no2= LS_Crack->no2;
   //  xtensor::xVector<>&  no3= LS_Crack->no3;
   xtensor::xTensor2<>& ri = LS_Crack->n123_row;
   xtensor::xTensor2<>& le = LS_Crack->n123_col;
   xtensor::xTensor2<>& rio = LS_Crack->no123_row;
   xtensor::xTensor2<>& leo = LS_Crack->no123_col;
   xtensor::xTensor2<>& curv0 = LS_Crack->curv1;
   xtensor::xTensor2<>& curv1 = LS_Crack->curv2;
   xtensor::xTensor2<>& curvo0 = LS_Crack->curvo1;
   xtensor::xTensor2<>& curvo1 = LS_Crack->curvo2;
   xtensor::xTensor2<>& curvo2 = LS_Crack->curvo3;
   xtensor::xVector<> dalpha = dalpha2;
// test : on ne prend que la composante radiale . . .
#if QRADIAL == 1
   for (int i = 0; i < 3; i++)
   {
      dalpha(i) = dalpha(i) - dalpha(i) * no3(i);
   }
#endif
   xtensor::xVector<> q(no1 * alpha);
   xtensor::xTensor2<> dq(curvo0);  // a vérifier!!
   dq *= alpha;
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         dq(i, j) += no1(i) * dalpha(j);
      }
   }

   // computation of Jhh
   xtensor::xTensor2<> GradDisph;
   grad_disp(geo_elem, geo_elem, GradDisph);

   //  GradDisph(2,2)=1;
   if (debug) cout << "Graddisph before transpose" << endl;
   if (debug) cout << GradDisph << endl;
   // Greg
   //  Transpo(GradDisph);
   if (debug) cout << "Graddisph after transpose" << endl;
   if (debug) cout << GradDisph << endl;
   xtensor::xTensor2<> Stressh = Phys * GradDisph;
   if (debug) cout << "Stressh" << endl;
   if (debug) cout << Stressh << endl;
   xtensor::xTensor2<> ShjiUhim(Stressh * GradDisph);

   if (debug) cout << "q" << endl;
   if (debug) cout << q << endl;

   if (debug) cout << "dq" << endl;
   if (debug) cout << dq << endl;

   double engh = 0.5 * ShjiUhim.trace();
   if (debug) cout << "engh " << engh << endl;
   xtensor::xTensor2<> eshelbyh;
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         eshelbyh(i, j) = -ShjiUhim(j, i);
         if (i == j) eshelbyh(i, i) += engh;
      }
   }

   Jhh = -Weight * DetJac * xfem::xFormProducts::product(dq, eshelbyh);
   vol = Weight * DetJac;

   // for the auxiliary modes
   static xtensor::xTensor2<> Stressl, Epsl, Gradl0Stress, Gradl1Stress;
   static xtensor::xVector<> Displ, Gradl0Disp, Gradl1Disp, Gradl00Disp, Gradl01Disp, Gradl11Disp;
   static xtensor::xTensor2<> Grad0Stress, Grad1Stress, Grad2Stress, GradDisp;
   static xtensor::xTensor2<> Grad0GradDisp, Grad1GradDisp, Grad2GradDisp;
   static xtensor::xTensor2<> temp0, temp1, temp2;
   static xtensor::xTensor2<> eshelbyIh;
   static double workIh;

   for (int mode = 0; mode < 3; mode++)
   {
      // the little l indicates that all these values are in the
      // local (level set) coordinate system

      GetMech3DAuxFields(LS_Crack->local, MatInfo, Stressl, Epsl, Displ, Gradl0Disp, Gradl1Disp, Gradl0Stress, Gradl1Stress,
                         Gradl00Disp, Gradl11Disp, Gradl01Disp, mode + 1);

      if (debug)
      {
         cout << "mode " << mode << " Stressl" << endl << Stressl << endl;
         cout << "mode " << mode << "Epsl" << endl << Epsl << endl;
         cout << "Displ" << endl << Displ << endl;
         cout << "Gradl0Disp" << endl << Gradl0Disp << endl;
         cout << "Gradl1Disp" << endl << Gradl1Disp << endl;
         cout << "Gradl0Stress" << endl << Gradl0Stress << endl;
         cout << "Gradl1Stress" << endl << Gradl1Stress << endl;
         cout << "Gradl00Disp" << endl << Gradl00Disp << endl;
         cout << "Gradl11Disp" << endl << Gradl11Disp << endl;
         cout << "Gradl01Disp" << endl << Gradl01Disp << endl;
         cout << " mode " << mode << " divStress is " << Gradl0Stress(0, 0) + Gradl1Stress(0, 1) << " "
              << Gradl0Stress(1, 0) + Gradl1Stress(1, 1) << " " << Gradl0Stress(2, 0) + Gradl1Stress(2, 1) << endl;
         cout << " mode " << mode << " must be zero " << 2. * Gradl00Disp(0) + Gradl11Disp(0) + Gradl01Disp(1) << " "
              << 2. * Gradl11Disp(1) + Gradl00Disp(1) + Gradl01Disp(0) << " " << Gradl00Disp(2) + Gradl11Disp(2) << endl;
      }

      xtensor::xTensor2<> Eps(leo * (Epsl * rio));
      xtensor::xTensor2<> Stress(leo * (Stressl * rio));

      if (debug)
      {
         cout << "Stress " << endl << Stress << endl;
         cout << "Eps" << endl << Eps << endl;
      }

      // we now change the local derivatives into global derivatives
      // page 38 BookI
      for (int a = 0; a < 3; a++)
      {
         for (int b = 0; b < 3; b++)
         {
            Grad0Stress(a, b) = Gradl0Stress(a, b) * ri(0, 0) + Gradl1Stress(a, b) * ri(1, 0);
            Grad1Stress(a, b) = Gradl0Stress(a, b) * ri(0, 1) + Gradl1Stress(a, b) * ri(1, 1);
            Grad2Stress(a, b) = Gradl0Stress(a, b) * ri(0, 2) + Gradl1Stress(a, b) * ri(1, 2);
         }
      }

      for (int a = 0; a < 3; a++)
      {
         GradDisp(a, 0) = Gradl0Disp(a) * ri(0, 0) + Gradl1Disp(a) * ri(1, 0);
         GradDisp(a, 1) = Gradl0Disp(a) * ri(0, 1) + Gradl1Disp(a) * ri(1, 1);
         GradDisp(a, 2) = Gradl0Disp(a) * ri(0, 2) + Gradl1Disp(a) * ri(1, 2);
      }
      if (debug) cout << "GradDisp" << endl << GradDisp << endl;

      if (debug)
      {
         cout << "mode is " << mode << endl;
         cout << "checking the symmetry of Grad0Stress " << endl << Grad0Stress << endl;
         cout << "checking the symmetry of Grad1Stress " << endl << Grad1Stress << endl;
         cout << "checking the symmetry of Grad2Stress " << endl << Grad2Stress << endl;
      }

      // xtensor::xTensor2<> curv0sym = curv0; curv0sym.symmetrize();
      // xtensor::xTensor2<> curv1sym = curv1; curv1sym.symmetrize();
      for (int a = 0; a < 3; a++)
      {
         for (int i = 0; i < 3; i++)
         {
            Grad0GradDisp(a, i) = le(i, 0) * Gradl00Disp(a) * ri(0, 0) + le(i, 1) * Gradl11Disp(a) * ri(1, 0) +
                                  le(i, 0) * Gradl01Disp(a) * ri(1, 0) + le(i, 1) * Gradl01Disp(a) * ri(0, 0) +
                                  Gradl0Disp(a) * curv0(i, 0) + Gradl1Disp(a) * curv1(i, 0);
            Grad1GradDisp(a, i) = le(i, 0) * Gradl00Disp(a) * ri(0, 1) + le(i, 1) * Gradl11Disp(a) * ri(1, 1) +
                                  le(i, 0) * Gradl01Disp(a) * ri(1, 1) + le(i, 1) * Gradl01Disp(a) * ri(0, 1) +
                                  Gradl0Disp(a) * curv0(i, 1) + Gradl1Disp(a) * curv1(i, 1);
            Grad2GradDisp(a, i) = le(i, 0) * Gradl00Disp(a) * ri(0, 2) + le(i, 1) * Gradl11Disp(a) * ri(1, 2) +
                                  le(i, 0) * Gradl01Disp(a) * ri(1, 2) + le(i, 1) * Gradl01Disp(a) * ri(0, 2) +
                                  Gradl0Disp(a) * curv0(i, 2) + Gradl1Disp(a) * curv1(i, 2);
         }
      }

      if (debug)
      {
         cout << "mode is " << mode << endl;
         for (int i = 0; i < 3; i++)
         {
            cout << "checking the symmetry of GradGradDisp i = " << i << endl
                 << Grad0GradDisp(i, 0) << " " << Grad0GradDisp(i, 1) << " " << Grad0GradDisp(i, 2) << endl
                 << Grad1GradDisp(i, 0) << " " << Grad1GradDisp(i, 1) << " " << Grad1GradDisp(i, 2) << endl
                 << Grad2GradDisp(i, 0) << " " << Grad2GradDisp(i, 1) << " " << Grad2GradDisp(i, 2) << endl;
         }
      }

      // we now transform the classical derivatives into
      // covariant derivatives
      // page 36 and 37 BookI
      xtensor::xTensor2<> CGradDisp(leo * GradDisp);
      if (debug) cout << " GradDisp avant correction \n" << GradDisp << endl;
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            CGradDisp(i, j) += Displ(0) * curvo0(i, j) + Displ(1) * curvo1(i, j) + Displ(2) * curvo2(i, j);
         }
      }
      if (debug)
      {
         cout << " leo    " << endl << leo << endl;
         cout << " curvo0 " << endl << curvo0 << endl;
         cout << " curvo1 " << endl << curvo1 << endl;
         cout << " curvo2 " << endl << curvo2 << endl;
         cout << " GradDisp " << GradDisp << endl;
         cout << " CGradDisp " << CGradDisp << endl;
      }
      ///////////////////////
      xtensor::xTensor2<> CGrad0CGradDisp(leo * Grad0GradDisp);
      xtensor::xTensor2<> CGrad1CGradDisp(leo * Grad1GradDisp);
      xtensor::xTensor2<> CGrad2CGradDisp(leo * Grad2GradDisp);
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            CGrad0CGradDisp(i, j) += curvo0(i, j) * GradDisp(0, 0) + curvo1(i, j) * GradDisp(1, 0) +
                                     curvo2(i, j) * GradDisp(2, 0) + curvo0(i, 0) * GradDisp(0, j) +
                                     curvo1(i, 0) * GradDisp(1, j) + curvo2(i, 0) * GradDisp(2, j);
            CGrad1CGradDisp(i, j) += curvo0(i, j) * GradDisp(0, 1) + curvo1(i, j) * GradDisp(1, 1) +
                                     curvo2(i, j) * GradDisp(2, 1) + curvo0(i, 1) * GradDisp(0, j) +
                                     curvo1(i, 1) * GradDisp(1, j) + curvo2(i, 1) * GradDisp(2, j);
            CGrad2CGradDisp(i, j) += curvo0(i, j) * GradDisp(0, 2) + curvo1(i, j) * GradDisp(1, 2) +
                                     curvo2(i, j) * GradDisp(2, 2) + curvo0(i, 2) * GradDisp(0, j) +
                                     curvo1(i, 2) * GradDisp(1, j) + curvo2(i, 2) * GradDisp(2, j);
         }
      }

      if (debug)
      {
         cout << "mode is " << mode << endl;
         for (int i = 0; i < 3; i++)
         {
            cout << "checking the symmetry of CGradCGradDisp i = " << i << endl
                 << CGrad0CGradDisp(i, 0) << " " << CGrad0CGradDisp(i, 1) << " " << CGrad0CGradDisp(i, 2) << endl
                 << CGrad1CGradDisp(i, 0) << " " << CGrad1CGradDisp(i, 1) << " " << CGrad1CGradDisp(i, 2) << endl
                 << CGrad2CGradDisp(i, 0) << " " << CGrad2CGradDisp(i, 1) << " " << CGrad2CGradDisp(i, 2) << endl;
         }
      }

      for (int i = 0; i < 3; i++)
      {
         temp0(i, 0) = curvo0(i, 0);
         temp0(i, 1) = curvo1(i, 0);
         temp0(i, 2) = curvo2(i, 0);
         temp1(i, 0) = curvo0(i, 1);
         temp1(i, 1) = curvo1(i, 1);
         temp1(i, 2) = curvo2(i, 1);
         temp2(i, 0) = curvo0(i, 2);
         temp2(i, 1) = curvo1(i, 2);
         temp2(i, 2) = curvo2(i, 2);
      }

      xtensor::xTensor2<> Stressrio(Stress * rio);
      xtensor::xTensor2<> leoStress(leo * Stress);
      xtensor::xTensor2<> CGrad0Stress(leo * (Grad0Stress * rio));
      CGrad0Stress += temp0 * Stressrio;
      CGrad0Stress += leoStress * Transpo(temp0);
      xtensor::xTensor2<> CGrad1Stress(leo * (Grad1Stress * rio));
      CGrad1Stress += temp1 * Stressrio;
      CGrad1Stress += leoStress * Transpo(temp1);
      xtensor::xTensor2<> CGrad2Stress(leo * (Grad2Stress * rio));
      CGrad2Stress += temp2 * Stressrio;
      CGrad2Stress += leoStress * Transpo(temp2);

      if (debug)
      {
         cout << "mode is " << mode << endl;
         cout << "checking the symmetry of CGrad0Stress " << endl << CGrad0Stress << endl;
         cout << "checking the symmetry of CGrad1Stress " << endl << CGrad1Stress << endl;
         cout << "checking the symmetry of CGrad2Stress " << endl << CGrad2Stress << endl;
      }

      // now eshelby interaction tensor
      // page 32 BookI
      // Stressh.contract(Eps) delta  - transpose ( Stressh * CGradDisp)
      //                             - transpose ( Stressl * GradDisph);

      ///////////////////
      //
      // finite element interaction integral
      //

      eshelbyIh = (Stressh * CGradDisp);
      if (debug) cout << " stressh gradisp " << eshelbyIh << endl;
      eshelbyIh += (Stress * GradDisph);
      if (debug) cout << " stress  gradisph " << (Stress * GradDisph) << endl;
      Transpo(eshelbyIh);
      eshelbyIh *= -1.0;
      workIh = xfem::xFormProducts::product(Stressh, Eps);
      if (debug) cout << " workIh " << workIh << endl;
      for (int i = 0; i < 3; i++) eshelbyIh(i, i) += workIh;

      xtensor::xVector<> div_eshelbyIh;
      div_eshelbyIh(0) = xfem::xFormProducts::product(CGrad0Stress, GradDisph);  //- Stressh.contract(CGrad0CGradDisp);
      div_eshelbyIh(1) = xfem::xFormProducts::product(CGrad1Stress, GradDisph);  //- Stressh.contract(CGrad1CGradDisp);
      div_eshelbyIh(2) = xfem::xFormProducts::product(CGrad2Stress, GradDisph);  //- Stressh.(CGrad2CGradDisp);
      for (int m = 0; m < 3; m++)
      {
         for (int i = 0; i < 3; i++)
         {
            div_eshelbyIh(m) -= Stressh(i, 0) * CGrad0CGradDisp(i, m) + Stressh(i, 1) * CGrad1CGradDisp(i, m) +
                                Stressh(i, 2) * CGrad2CGradDisp(i, m);
         }
      }
      xtensor::xVector<> divStress;
      for (int i = 0; i < 3; i++)
      {
         divStress(i) = CGrad0Stress(i, 0) + CGrad1Stress(i, 1) + CGrad2Stress(i, 2);
      }
      div_eshelbyIh -= (divStress * GradDisph);

      // cout << "mode " << mode << "divStress"  << divStress << endl;
      // cout << "div_eshelbyIh"  << div_eshelbyIh << endl;

      // Ih(mode) = - Weight * DetJac * (dq.contract(eshelbyIh) + (q * div_eshelbyIh) );
      // not corrected works best so far ...
      // orig Ih(mode) = - Weight * DetJac * dq.contract(eshelbyIh);
      if (debug)
         cout << "mode " << mode << " xfem::xFormProducts::product(dq, eshelbyIh); "
              << xfem::xFormProducts::product(dq, eshelbyIh) << endl;
      Ih(mode) = -Weight * DetJac * xfem::xFormProducts::product(dq, eshelbyIh);
//    Ih(mode) = - Weight * DetJac *  xfem::xForm::product(dq, eshelbyIh) + (q * div_eshelbyIh) ;
#if 1  // no computation of the mode mode interaction
       // to speed up
      //
      // mode mode interaction
      //
      xtensor::xTensor2<> eshelbyI(Stress * CGradDisp);
      double workI = eshelbyI.trace();
      Transpo(eshelbyI);
      eshelbyI *= -1.0;
      for (int i = 0; i < 3; i++) eshelbyI(i, i) += 0.5 * workI;

      // see page 39 bookI
      // in fact we divide we divide by two the results of page 39
      xtensor::xVector<> div_eshelbyI;
      div_eshelbyI(0) = xfem::xFormProducts::product(CGrad0Stress, CGradDisp);  //- Stress.contract(CGrad0CGradDisp);
      div_eshelbyI(1) = xfem::xFormProducts::product(CGrad1Stress, CGradDisp);  //- Stress.contract(CGrad1CGradDisp);
      div_eshelbyI(2) = xfem::xFormProducts::product(CGrad2Stress, CGradDisp);  //- Stress.contract(CGrad2CGradDisp);
      for (int m = 0; m < 3; m++)
      {
         for (int i = 0; i < 3; i++)
         {
            div_eshelbyI(m) -= Stress(i, 0) * CGrad0CGradDisp(i, m) + Stress(i, 1) * CGrad1CGradDisp(i, m) +
                               Stress(i, 2) * CGrad2CGradDisp(i, m);
         }
      }
      div_eshelbyI -= (divStress * CGradDisp);

      // cout << "mode " << mode << "div_eshelbyI"  << div_eshelbyI << endl;

      // I(mode, mode) = - Weight * DetJac * (dq.contract(eshelbyI) + (q * div_eshelbyI) );
      // not corrected works best so far ...
      I(mode, mode) = -Weight * DetJac * xfem::xFormProducts::product(dq, eshelbyI);

#endif
   }
   if (debug)
   {
      cout << "Jhh_l " << Jhh << endl;
      cout << "Ih_l  " << Ih << endl;
      cout << "I_l   " << endl;
      cout << I << endl;
      cout << "vol_l " << vol << endl;
   }
   return;
}

void JGrid3D(const xtensor::xPoint& pointJfront, double rx, double ry, double rz, const xtensor::xTensor2<>& mat,
             std::map<int, xtensor::xPoint>& JCoord, std::map<int, std::vector<int>>& JElements)
{
   // Nb of cells in the J-Domain Box
   const int CELLX = 2;
   const int CELLY = 2;
   const int CELLZ = 2;
   //
   int cellx = CELLX, celly = CELLY, cellz = CELLZ;
   double dx = rx / (double)cellx, dy = ry / (double)celly, dz = rz / (double)cellz;
   //   int    celltotal = cellx*celly*cellz;
   int ndx = cellx + 1, ndy = celly + 1, ndz = cellz + 1;
   //   int    ndtotal = ndx*ndy*ndz;
   int i, j, k;

   xtensor::xPoint xyz, xyzrot;
   // typedef std::map<int,double[3]>::value_type coord_add;

   // clear the maps
   JCoord.clear();
   JElements.clear();

#if DEBUGSUKU >= 1
   // print the inverse (transpose) transformation matrix
   for (i = 0; i < 3; ++i)
   {
      printf("Row %d of mat is: %22.15e %22.15e %22.15e \n", i, mat(i, 0), mat(i, 1), mat(i, 2));
   }
   printf("J evaluated at the point: %22.15e %22.15e %22.15e \n", pointJfront[0], pointJfront[1], pointJfront[2]);
#endif

   int id = 1;
   for (k = 0; k < ndz; ++k)
   {
      for (j = 0; j < ndy; ++j)
      {
         for (i = 0; i < ndx; ++i)
         {
            xyz(0) = i * dx - 0.5 * rx;
            xyz(1) = j * dy - 0.5 * ry;
            xyz(2) = k * dz - 0.5 * rz;
            // rotate the coordinates using mat[i][j]
            for (int i1 = 0; i1 < 3; ++i1)
            {
               xyzrot(i1) = 0.0;
               for (int j1 = 0; j1 < 3; ++j1)
               {
                  xyzrot(i1) += mat(i1, j1) * xyz(j1);
               }
               // translation
               xyzrot(i1) += pointJfront(i1);
            }
            // xyzrot is the global coordinate
            JCoord.insert(make_pair(id, xyzrot));
#if DEBUGSUKU >= 1
            printf("JDomain Grid Info\n");
            printf("Local  Coord %d %d %d is: %22.15e %22.15e %22.15e \n", i, j, k, xyz[0], xyz[1], xyz[2]);
            printf("Global Coord %d %d %d is: %22.15e %22.15e %22.15e\n ", i, j, k, xyzrot[0], xyzrot[1], xyzrot[2]);
#endif
            ++id;
         }
      }
   }

   //
   // generate cell connectivities
   //
   // typedef std::map<int,int[8]>::value_type connect_add;
   std::vector<int> connect(8);

   int index, k1;
   int nadd = ndx * ndy;

   int node_id = 1;
   id = 1;
   for (k = 0; k < cellz; ++k)
   {
      for (j = 0; j < celly; ++j)
      {
         index = j * ndx + k * nadd;
         for (i = 0; i < cellx; ++i)
         {
            k1 = index + i;
            connect[0] = k1 + node_id;
            connect[1] = connect[0] + 1;
            connect[2] = connect[1] + ndx;
            connect[3] = connect[2] - 1;
            connect[4] = connect[0] + nadd;
            connect[5] = connect[1] + nadd;
            connect[6] = connect[2] + nadd;
            connect[7] = connect[3] + nadd;
            JElements.insert(make_pair(id, connect));
#if DEBUGSUKU >= 1
            printf("Connectivity for element %d :", id);
            for (int i1 = 0; i1 < 8; ++i1)
            {
               printf("%d ", connect[i1]);
            }
            printf("\n");
#endif
            ++id;
         }
      }
   }
}

void GetGradqLoc(int postid, const xtensor::xPoint& xyz, xtensor::xTensor2<>& GradqLoc, xtensor::xPoint& xyzloc,
                 const std::map<int, xtensor::xVector<>>& JBoxDimension, const std::map<int, xtensor::xPoint>& JFrontPoint,
                 const std::map<int, xtensor::xTensor2<>>& JMap)
{
   // compute gradient of q in the local crack front coordinate system
   std::map<int, xtensor::xVector<>>::const_iterator it3End = JBoxDimension.end();
   std::map<int, xtensor::xVector<>>::const_iterator itfind3 = JBoxDimension.find(postid);
   if (itfind3 == it3End)
   {
      assert(1 == 0);
   }

   xtensor::xVector<> boxdim = itfind3->second;
   xtensor::xVector<> scale = boxdim * 0.5;

   std::map<int, xtensor::xPoint>::const_iterator itf3End = JFrontPoint.end();
   std::map<int, xtensor::xPoint>::const_iterator itffind3 = JFrontPoint.find(postid);
   if (itffind3 == itf3End)
   {
      assert(1 == 0);
   }

   xtensor::xPoint frontpoint = itffind3->second;

   std::map<int, xtensor::xTensor2<>>::const_iterator itmEnd = JMap.end();
   std::map<int, xtensor::xTensor2<>>::const_iterator itmf = JMap.find(postid);
   if (itmf == itmEnd)
   {
      assert(1 == 0);
   }

   xtensor::xTensor2<> rotmat = itmf->second;

   int i, j;
   for (i = 0; i < 3; ++i)
   {
      for (j = 0; j < 3; ++j)
      {
         GradqLoc(i, j) = 0.0;
      }
   }

   xtensor::xVector<> pointglo(frontpoint, xyz);
   //  double pointglo[3] = { xyz[0] - frontpoint[0], xyz[1] - frontpoint[1], xyz[2] - frontpoint[2] };
   // xyzloc = rotmat * pointglo;
   for (i = 0; i < 3; ++i)
   {
      xyzloc(i) = 0.;
      for (j = 0; j < 3; ++j)
      {
         xyzloc(i) += rotmat(i, j) * pointglo(j);
      }
   }
#if DEBUGNIC >= 1
   printvecwithname("xyz", xyz);
   printvecwithname("xyzloc", xyzloc);
   printvecwithname("frontpoint", frontpoint);
   printmatwithname("Inside Gradqloc rotmat is ", rotmat);
#endif

   xtensor::xPoint xyzprime;
   for (i = 0; i < 3; ++i)
   {
      xyzprime(i) = fabs(xyzloc(i) / scale(i));
   }

   if (xyzprime(0) <= 1.0 && xyzprime(1) <= 1.0 && xyzprime(2) <= 1.0)
   {
#if QMODIF == 3
      if (xyzprime(0) <= TCUBE && xyzprime(1) <= TCUBE)
      {
         // gradient of q in the local coordinate system
         GradqLoc(0, 0) = 0.0;
         GradqLoc(0, 1) = 0.0;
         GradqLoc(0, 2) = 0.0;
      }
      else if (xyzprime(0) <= TCUBE)
      {
         // gradient of q in the local coordinate system
         GradqLoc(0, 0) = 0.0;
         GradqLoc(0, 1) = -1.0 * ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / (scale(1) * (1 - TCUBE));
         GradqLoc(0, 2) = 0.0;
      }
      else if (xyzprime(1) <= TCUBE)
      {
         // gradient of q in the local coordinate system
         GradqLoc(0, 0) = -1.0 * ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / (scale(0) * (1 - TCUBE));
         GradqLoc(0, 1) = 0.0;
         GradqLoc(0, 2) = 0.0;
      }
      else
      {
         GradqLoc(0, 0) = -((1. - xyzprime(1)) / (1 - TCUBE)) * ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / (scale(0) * (1 - TCUBE));
         GradqLoc(0, 1) = -((1. - xyzprime(0)) / (1 - TCUBE)) * ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / (scale(1) * (1 - TCUBE));
         GradqLoc(0, 2) = -0.0;
         ;
      }

#elif QMODIF == 2
      if (xyzprime(0) <= TCUBE && xyzprime(1) <= TCUBE && xyzprime(2) <= TCUBE)
      {
         GradqLoc(0, 0) = 0.0;
         GradqLoc(0, 1) = 0.0;
         GradqLoc(0, 2) = 0.0;
      }
      else if (xyzprime(0) <= TCUBE && xyzprime(1) <= TCUBE)
      {
         GradqLoc(0, 0) = 0.0;
         GradqLoc(0, 1) = 0.0;
         GradqLoc(0, 2) = -1 * ((xyzloc(2) >= 0.0) ? 1.0 : -1.0) / (scale(2) * (1 - TCUBE));
      }
      else if (xyzprime(0) <= TCUBE && xyzprime(2) <= TCUBE)
      {
         GradqLoc(0, 0) = 0.0;
         GradqLoc(0, 1) = -1 * ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / (scale(1) * (1 - TCUBE));
         GradqLoc(0, 2) = 0.0;
      }
      else if (xyzprime(1) <= TCUBE && xyzprime(2) <= TCUBE)
      {
         GradqLoc(0, 0) = -1 * ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / (scale(0) * (1 - TCUBE));
         GradqLoc(0, 1) = 0.0;
         GradqLoc(0, 2) = 0.0;
      }
      else if (xyzprime(0) <= TCUBE)
      {
         GradqLoc(0, 0) = 0.0;
         GradqLoc(0, 1) = -((1. - xyzprime(2)) / (1 - TCUBE)) * ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / (scale(1) * (1 - TCUBE));
         GradqLoc(0, 2) = -((1. - xyzprime(1)) / (1 - TCUBE)) * ((xyzloc(2) >= 0.0) ? 1.0 : -1.0) / (scale(2) * (1 - TCUBE));
      }
      else if (xyzprime(1) <= TCUBE)
      {
         GradqLoc(0, 0) = -((1. - xyzprime(2)) / (1 - TCUBE)) * ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / (scale(0) * (1 - TCUBE));
         GradqLoc(0, 1) = 0.0;
         GradqLoc(0, 2) = -((1. - xyzprime(0)) / (1 - TCUBE)) * ((xyzloc(2) >= 0.0) ? 1.0 : -1.0) / (scale(2) * (1 - TCUBE));
      }
      else if (xyzprime(2) <= TCUBE)
      {
         GradqLoc(0, 0) = -((1. - xyzprime(1)) / (1 - TCUBE)) * ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / (scale(0) * (1 - TCUBE));
         GradqLoc(0, 1) = -((1. - xyzprime(0)) / (1 - TCUBE)) * ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / (scale(1) * (1 - TCUBE));
         GradqLoc(0, 2) = 0.0;
      }
      else
      {
         GradqLoc(0, 0) = -((1. - xyzprime(1)) / (1 - TCUBE)) * ((1. - xyzprime(2)) / (1 - TCUBE)) *
                          ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / (scale(0) * (1 - TCUBE));
         GradqLoc(0, 1) = -((1. - xyzprime(0)) / (1 - TCUBE)) * ((1. - xyzprime(2)) / (1 - TCUBE)) *
                          ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / (scale(1) * (1 - TCUBE));
         GradqLoc(0, 2) = -((1. - xyzprime(0)) / (1 - TCUBE)) * ((1. - xyzprime(1)) / (1 - TCUBE)) *
                          ((xyzloc(2) >= 0.0) ? 1.0 : -1.0) / (scale(2) * (1 - TCUBE));
      }
#elif QMODIF == 1
      double xz = (1.0 - xyzprime(0));
      double yz = (1.0 - xyzprime(1));
      // double qw = (1.0 - xyzprime(0))*(1.0 - xyzprime(1)); -Wunused-but-set-variable
      // gradient of q in the local coordinate system
      GradqLoc(0, 0) = -yz * ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / scale(0);
      GradqLoc(0, 1) = -xz * ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / scale(1);
      GradqLoc(0, 2) = 0.0;
#else
      double xy = (1.0 - xyzprime(0)) * (1.0 - xyzprime(1));
      double xz = (1.0 - xyzprime(0)) * (1.0 - xyzprime(2));
      double yz = (1.0 - xyzprime(1)) * (1.0 - xyzprime(2));
      // double qw = (1.0 - xyzprime(0))*(1.0 - xyzprime(1))*(1. - xyzprime(2)); -Wunused-but-set-variable
      // gradient of q in the local coordinate system
      GradqLoc(0, 0) = -yz * ((xyzloc(0) >= 0.0) ? 1.0 : -1.0) / scale(0);
      GradqLoc(0, 1) = -xz * ((xyzloc(1) >= 0.0) ? 1.0 : -1.0) / scale(1);

      GradqLoc(0, 2) = -xy * ((xyzloc(2) >= 0.0) ? 1.0 : -1.0) / scale(2);
#endif
   }

   else
   {
      //    cerr << "Integration point" << xyz(0) << " " << xyz(1) << " " << xyz(2) <<
      //      " in JIntegral calculation is outside JDomain BOX \n";
      //    cerr << "xyzprime" << xyzprime(0) << " " << xyzprime(1) << " " << xyzprime(2) << endl;
      assert(1 == 0);
   }

#if DEBUGSUKU >= 1
   printf("Point on the crack front is: %22.15e %22.15e %22.15e \n", frontpoint(0), frontpoint(1), frontpoint(2));
   printf("Integration point is   : %22.15e %22.15e %22.15e \n", xyz(0), xyz(1), xyz(2));
   printf("Local position is: %22.15e %22.15e %22.15e \n", xyzloc(0), xyzloc(1), xyzloc(2));
   printmatwithname("rotmat", rotmat);
   printf("GradqLoc row 1 is: %22.15e %22.15e and %22.15e \n", GradqLoc(0, 0), GradqLoc(0, 1), GradqLoc(0, 2));
   printf("GradqLoc row 2 is: %22.15e %22.15e and %22.15e \n", GradqLoc(1, 0), GradqLoc(1, 1), GradqLoc(1, 2));
   printf("GradqLoc row 3 is: %22.15e %22.15e and %22.15e \n", GradqLoc(2, 0), GradqLoc(2, 1), GradqLoc(2, 2));
#endif

   return;
}

double AngleWithMaxHoopStress(const double& K1, const double& K2)
{
   const bool debug = true;
   double rat, thetac1, thetac2, hoop1;
   rat = K1 / K2;
   thetac1 = 2. * atan(0.25 * (rat + sqrt(rat * rat + 8.)));
   thetac2 = 2. * atan(0.25 * (rat - sqrt(rat * rat + 8.)));
   //
   //  printf("theta1 is %12.5e and theta2 is %12.5e\n", thetac1, thetac2);
   //
   hoop1 = K1 * (3. * cos(thetac1 / 2.) + cos(3. * thetac1 / 2.)) - 3. * K2 * (sin(thetac1 / 2.) + sin(3. * thetac1 / 2.));
   //    - K2*(3.*sin(thetac1/2.) + sin(3.*thetac1/2.));//warning
   double hoop2 = K1 * (3. * cos(thetac2 / 2.) + cos(3. * thetac2 / 2.)) - 3. * K2 * (sin(thetac2 / 2.) + sin(3. * thetac2 / 2.));
   if (debug)
   {
      cout << " K1 = " << K1 << endl;
      cout << " K2 = " << K2 << endl;
      cout << " rat (K1/K2) = " << rat << endl;
      cout << " thetac1 = " << thetac1 << endl;
      cout << " thetac2 = " << thetac2 << endl;
      cout << " hoop1 : " << hoop1 << endl;
      cout << " hoop2 : " << hoop2 << endl;
   }
   if ((hoop1) >= (hoop2))
      return thetac1;
   else
      return thetac2;
}

double GetqLoc(int postid, const xtensor::xPoint& xyz, const std::map<int, xtensor::xVector<>>& JBoxDimension,
               const std::map<int, xtensor::xPoint>& JFrontPoint, const std::map<int, xtensor::xTensor2<>>& JMap)
{
   // compute q in the local crack front coordinate system

   std::map<int, xtensor::xVector<>>::const_iterator it3End = JBoxDimension.end();
   std::map<int, xtensor::xVector<>>::const_iterator itfind3 = JBoxDimension.find(postid);
   if (itfind3 == it3End)
   {
      assert(1 == 0);
   }

   xtensor::xVector<> boxdim = itfind3->second;

   xtensor::xVector<> scale = boxdim * 0.5;
   //   double scale(3];
   //   scale(0] = 0.5*boxdim(0];
   //   scale(1] = 0.5*boxdim(1];
   //   scale(2] = 0.5*boxdim(2];

   std::map<int, xtensor::xPoint>::const_iterator itf3End = JFrontPoint.end();
   std::map<int, xtensor::xPoint>::const_iterator itffind3 = JFrontPoint.find(postid);
   if (itffind3 == itf3End)
   {
      assert(1 == 0);
   }

   xtensor::xPoint frontpoint = itffind3->second;

   std::map<int, xtensor::xTensor2<>>::const_iterator itmEnd = JMap.end();
   std::map<int, xtensor::xTensor2<>>::const_iterator itmf = JMap.find(postid);
   if (itmf == itmEnd)
   {
      assert(1 == 0);
   }

   xtensor::xTensor2<> rotmat = itmf->second;

   int i;

   xtensor::xVector<> pointloc(frontpoint, xyz);
   xtensor::xVector<> xyzloc = rotmat * pointloc;

   // double pointloc[3] = { xyz[0] - frontpoint[0], xyz[1] - frontpoint[1], xyz[2] - frontpoint[2] };
   //   double xyzloc[3];
   //   for (i = 0; i < 3; ++i) {
   //     xyzloc[i] = 0.;
   //     for (j = 0; j < 3; ++j) {
   //       xyzloc[i] += rotmat(i,j) * pointloc[j];
   //     }
   //   }

   xtensor::xPoint xyzprime;
   for (i = 0; i < 3; ++i)
   {
      xyzprime(i) = fabs(xyzloc(i) / scale(i));
   }

   double qw = 0.0;

   if (xyzprime(0) <= 1.0 && xyzprime(1) <= 1.0 && xyzprime(2) <= 1.0)
   {
#if QMODIF == 3
      if (xyzprime(0) <= TCUBE && xyzprime(1) <= TCUBE)
      {
         qw = 1.0;
      }
      else if (xyzprime(0) <= TCUBE)
      {
         qw = (1. - xyzprime(1)) / (1 - TCUBE);
      }
      else if (xyzprime(1) <= TCUBE)
      {
         qw = (1. - xyzprime(0)) / (1 - TCUBE);
      }
      else
      {
         qw = ((1. - xyzprime(0)) / (1 - TCUBE)) * ((1. - xyzprime(1)) / (1 - TCUBE));
      }

#elif QMODIF == 2
      if (xyzprime(0) <= TCUBE && xyzprime(1) <= TCUBE && xyzprime(2) <= TCUBE)
      {
         qw = 1.0;
      }
      else if (xyzprime(0) <= TCUBE && xyzprime(1) <= TCUBE)
      {
         qw = (1. - xyzprime(2)) / (1 - TCUBE);
      }
      else if (xyzprime(0) <= TCUBE && xyzprime(2) <= TCUBE)
      {
         qw = (1. - xyzprime(1)) / (1 - TCUBE);
      }
      else if (xyzprime(1) <= TCUBE && xyzprime(2) <= TCUBE)
      {
         qw = (1. - xyzprime(0)) / (1 - TCUBE);
      }
      else if (xyzprime(0) <= TCUBE)
      {
         qw = ((1. - xyzprime(1)) / (1 - TCUBE)) * ((1. - xyzprime(2)) / (1 - TCUBE));
      }
      else if (xyzprime(1) <= TCUBE)
      {
         qw = ((1. - xyzprime(0)) / (1 - TCUBE)) * ((1. - xyzprime(2)) / (1 - TCUBE));
      }
      else if (xyzprime(2) <= TCUBE)
      {
         qw = ((1. - xyzprime(0)) / (1 - TCUBE)) * ((1. - xyzprime(1)) / (1 - TCUBE));
      }
      else
      {
         qw = ((1. - xyzprime(0)) / (1 - TCUBE)) * ((1. - xyzprime(1)) / (1 - TCUBE)) * ((1. - xyzprime(2)) / (1 - TCUBE));
      }

#elif QMODIF == 1
      qw = (1.0 - xyzprime(0)) * (1.0 - xyzprime(1));
#else
      qw = (1.0 - xyzprime(0)) * (1.0 - xyzprime(1)) * (1. - xyzprime(2));
#endif
   }

#if DEBUGSUKU >= 1
   printf("q is: %22.15e \n", qw);
#endif

   return qw;
}

void CreateMeshWithPostProElements(xMesh* mesh, const std::map<int, xtensor::xPoint>& new_coords,
                                   const std::map<int, std::vector<int>>& new_conn3D)
{
   std::map<int, xtensor::xPoint>::const_iterator itc, itcEnd;
   std::map<int, std::vector<int>>::const_iterator it3d, it3dEnd;

   itc = new_coords.begin();
   itcEnd = new_coords.end();
   for (; itc != itcEnd; ++itc)
   {
      xtensor::xPoint p(itc->second);
      mesh->getMesh().createVertex(itc->first, p(0), p(1), p(2), nullptr);
   }

   it3d = new_conn3D.begin();
   it3dEnd = new_conn3D.end();
   for (; it3d != it3dEnd; ++it3d)
   {
      std::vector<int> pv = it3d->second;
      mesh->getMesh().createHexWithVertices(pv[0], pv[1], pv[2], pv[3], pv[4], pv[5], pv[6], pv[7], nullptr);
   }

   return;
}

}  // end namespace xcrack
