/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#include "JintModalParks.h"

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "JintParks.h"
#include "lCrack.h"
#include "xAlgorithm.h"
#include "xDebug.h"
#include "xExportGmsh.h"
#include "xMaterial.h"
#include "xMaterialManager.h"
#include "xMaterialSensitivity.h"
#include "xMesh.h"
#include "xTensor2.h"
#include "xTensor4.h"
#include "xVector.h"
#include "xcFormLinearEnergyRelease.h"
#include "xcValueJs.h"
//#include "xLinearSystemSolverSuperLU.h"
#include "mVertex.h"
#include "xTensorOperations.h"
#include "xcFilters.h"

namespace xcrack
{
/*
void lCrack::getJint3DModalParks(xField<>& disp_l,
                           xIntegrationRule& integrator_vol,
                           xIntegrationRule& integrator_bnd,
                           const double& rho_cylinder,
                            int nb_layers_cylinder,
                           int nb_layers_core,
                           int nb_modes,
                           std::string filename)
{
xEvalGradField<xtool::xIdentity<xtensor::xTensor2<> > > grad_disp(disp_l);
xUniformMaterialSensitivity<xtensor::xTensor4<> > hooke("strain");

xcEvalJintTermElasticity eval_terms(grad_disp, hooke, *this );

getJint3DModalParks(disp_l, integrator_vol, integrator_bnd,
                    rho_cylinder, nb_layers_cylinder, nb_layers_core, nb_modes,filename,
                    eval_terms);
};

void lCrack::getJint3DModalParksV2(xField<>& disp_l,
                           xIntegrationRule& integrator_vol,
                           xIntegrationRule& integrator_bnd,
                           const double& rho_cylinder,
                            int nb_layers_cylinder,
                           int nb_layers_core,
                           int nb_modes,
                           std::string filename)
{
xEvalGradField<xtool::xIdentity<xtensor::xTensor2<> > > grad_disp(disp_l);
xUniformMaterialSensitivity<xtensor::xTensor4<> > hooke("strain");

xcEvalJintTermElasticity eval_terms(grad_disp, hooke, *this );

getJint3DModalParksV2( integrator_vol, integrator_bnd,
                    rho_cylinder, nb_layers_cylinder, nb_layers_core, nb_modes,filename,
                    eval_terms);
};

*/

void exportSifs(const std::string& filename, const xField<>& j_modal, xRegion& frontPart)
{
   xEvalField<xtool::xIdentity<double>> eval(j_modal);

   cout << "entering printSifs" << endl;

   std::map<double, vector<double>> res;
   xIter it = frontPart.begin(1);
   xIter end = frontPart.end(1);
   xUpperCreatorRecursive integ2appro(3);
   for (; it != end; ++it)
   {
      mEntity* e_integ = *it;
      xfem::xGeomElem geom_integ(e_integ);
      mEntity* e_appro = integ2appro(e_integ);
      xfem::xGeomElem geom_appro(e_appro);
      for (int ipt = 0; ipt < 2; ++ipt)
      {
         double u = (ipt == 0) ? -1. : 1.;
         geom_integ.setUVW({u, 0, 0});
         if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
         xtensor::xPoint xyz = geom_integ.getXYZ();
         // double s = lss.getVal(e_appro, geom_appro.getUVW());
         double s = 0.;
         // lss_1D.getVal(e_integ, geom_integ.getUVW());

         double val_j, val_k1, val_k2, val_k3;
         xcValueJsAndSifs::js(1);
         eval(&geom_appro, &geom_integ, val_j);
         xcValueJsAndSifs::sifs(1);
         eval(&geom_appro, &geom_integ, val_k1);
         xcValueJsAndSifs::sifs(2);
         eval(&geom_appro, &geom_integ, val_k2);
         xcValueJsAndSifs::sifs(3);
         eval(&geom_appro, &geom_integ, val_k3);

         std::vector<double> vals;
         vals.push_back(xyz(0));
         vals.push_back(xyz(1));
         vals.push_back(xyz(2));
         vals.push_back(val_k1);
         vals.push_back(val_k2);
         vals.push_back(val_k3);
         vals.push_back(val_j);

         res[s] = vals;
      }
   }

   string filename_s = filename + "_s.txt";
   std::ofstream out(filename_s.c_str());
   out.setf(ios_base::scientific);
   out.setf(ios_base::right);
   out.width(12);
   out << "s";
   out.width(12);
   out << "x";
   out.width(12);
   out << "y";
   out.width(12);
   out << "z";
   out.width(12);
   out << "K1";
   out.width(12);
   out << "K2";
   out.width(12);
   out << "K3";
   out.width(12);
   out << "J\n";

   std::map<double, vector<double>>::const_iterator its = res.begin();
   std::map<double, vector<double>>::const_iterator itse = res.end();
   for (; its != itse; ++its)
   {
      double s = its->first;
      vector<double> v = its->second;
      out.width(12);
      out << setprecision(3) << s;
      for (int i = 0; i < 7; ++i)
      {
         out.width(12);
         out << v[i];
      }
      out << endl;
      // copy(its->second.begin(), its->second.end(), std::ostream_iterator<double>(out, " "));
   }
   out.close();
   {
      string filename_modes = filename + "_modes.txt";
      std::ofstream out(filename_modes.c_str());
      out.setf(ios_base::scientific);
      out.setf(ios_base::right);

      out.setf(ios_base::scientific);
      out.setf(ios_base::right);
      out.width(12);
      out << "mode";
      out.width(12);
      out << "K1";
      out.width(12);
      out << "K2";
      out.width(12);
      out << "K3";
      out.width(12);
      out << "J\n";

      const xValueManagerDist<double>& double_manager = *j_modal.getValueManager();

      xIsPhys test_1d("J_modal");
      xValueManagerDist<double>::map_const_iterator itm = double_manager.begin();
      xValueManagerDist<double>::map_const_iterator itmEnd = double_manager.end();

      std::map<string, vector<double>> res;

      for (; itm != itmEnd; ++itm)
      {
         assert(test_1d(itm->first.getPhys()));
         string mode_name = xKeyInfo::getGeomName(itm->first.getGeom());

         xcValueJsAndSifs::sifs(1);
         double K1 = itm->second->getVal();
         xcValueJsAndSifs::sifs(2);
         double K2 = itm->second->getVal();
         xcValueJsAndSifs::sifs(3);
         double K3 = itm->second->getVal();
         xcValueJsAndSifs::js(1);
         double J = itm->second->getVal();

         std::vector<double> v;
         v.push_back(K1);
         v.push_back(K2);
         v.push_back(K3);
         v.push_back(J);
         res[mode_name] = v;
      }

      std::map<string, vector<double>>::const_iterator its = res.begin();
      std::map<string, vector<double>>::const_iterator itse = res.end();
      for (; its != itse; ++its)
      {
         string s = its->first;
         vector<double> v = its->second;
         out.width(12);
         out << setprecision(3) << s;
         for (size_t i = 0; i < v.size(); ++i)
         {
            out.width(12);
            out << v[i];
         }
         out << endl;
      }
      out.close();
   }
}

void lCrack::printSifsModal(const std::string& filename, const xField<>& j_modal, mEntity* seed)
{
   {
      xEvalField<xtool::xIdentity<double>> eval(j_modal);

      cout << "entering printSifs" << endl;

      std::map<double, vector<double>> res;
      xIter it = mesh_crack_front->begin(1);
      xIter end = mesh_crack_front->end(1);
      xUpperCreatorRecursive integ2appro(dim());
      if (seed == nullptr)
      {
         for (; it != end; ++it)
         {
            mEntity* e_integ = *it;
            xfem::xGeomElem geom_integ(e_integ);
            mEntity* e_appro = integ2appro(e_integ);
            xfem::xGeomElem geom_appro(e_appro);
            for (int ipt = 0; ipt < 2; ++ipt)
            {
               double u = (ipt == 0) ? -1. : 1.;
               geom_integ.setUVW({u, 0, 0});

               if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());

               xtensor::xPoint xyz = geom_integ.getXYZ();
               // double s = lss.getVal(e_appro, geom_appro.getUVW());
               double s = lss_1D.getVal(e_integ, geom_integ.getUVW());

               double val_j, val_k1, val_k2, val_k3;
               xcValueJsAndSifs::js(1);
               eval(&geom_appro, &geom_integ, val_j);
               xcValueJsAndSifs::sifs(1);
               eval(&geom_appro, &geom_integ, val_k1);
               xcValueJsAndSifs::sifs(2);
               eval(&geom_appro, &geom_integ, val_k2);
               xcValueJsAndSifs::sifs(3);
               eval(&geom_appro, &geom_integ, val_k3);

               std::vector<double> vals;
               vals.push_back(xyz(0));
               vals.push_back(xyz(1));
               vals.push_back(xyz(2));
               vals.push_back(val_k1);
               vals.push_back(val_k2);
               vals.push_back(val_k3);
               vals.push_back(val_j);

               res[s] = vals;
            }
         }
      }
      else
      {
         xtensor::xPoint xyz = static_cast<mVertex*>(seed)->point();
         double s = 0.;
         xfem::xGeomElem geom_integ(seed);
         mEntity* e_appro = integ2appro(seed);
         xfem::xGeomElem geom_appro(e_appro);
         geom_integ.setUVW({0, 0, 0});
         if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
         double val_j, val_k1, val_k2, val_k3;
         xcValueJsAndSifs::js(1);
         eval(&geom_appro, &geom_integ, val_j);
         xcValueJsAndSifs::sifs(1);
         eval(&geom_appro, &geom_integ, val_k1);
         xcValueJsAndSifs::sifs(2);
         eval(&geom_appro, &geom_integ, val_k2);
         xcValueJsAndSifs::sifs(3);
         eval(&geom_appro, &geom_integ, val_k3);

         std::vector<double> vals;
         vals.push_back(xyz(0));
         vals.push_back(xyz(1));
         vals.push_back(xyz(2));
         vals.push_back(val_k1);
         vals.push_back(val_k2);
         vals.push_back(val_k3);
         vals.push_back(val_j);
         res[s] = vals;
      }

      string filename_s = filename + "_s.txt";
      std::ofstream out(filename_s.c_str());
      out.setf(ios_base::scientific);
      out.setf(ios_base::right);
      out.width(12);
      out << "s";
      out.width(12);
      out << "x";
      out.width(12);
      out << "y";
      out.width(12);
      out << "z";
      out.width(12);
      out << "K1";
      out.width(12);
      out << "K2";
      out.width(12);
      out << "K3";
      out.width(12);
      out << "J\n";

      std::map<double, vector<double>>::const_iterator its = res.begin();
      std::map<double, vector<double>>::const_iterator itse = res.end();
      for (; its != itse; ++its)
      {
         double s = its->first;
         vector<double> v = its->second;
         out.width(12);
         out << setprecision(3) << s;
         for (int i = 0; i < 7; ++i)
         {
            out.width(12);
            out << v[i];
         }
         out << endl;
         // copy(its->second.begin(), its->second.end(), std::ostream_iterator<double>(out, " "));
      }
      out.close();
   }
   //

   {
      string filename_modes = filename + "_modes.txt";
      std::ofstream out(filename_modes.c_str());
      out.setf(ios_base::scientific);
      out.setf(ios_base::right);

      out.setf(ios_base::scientific);
      out.setf(ios_base::right);
      out.width(12);
      out << "mode";
      out.width(12);
      out << "K1";
      out.width(12);
      out << "K2";
      out.width(12);
      out << "K3";
      out.width(12);
      out << "J\n";

      const xValueManagerDist<double>& double_manager = *j_modal.getValueManager();

      xIsPhys test_1d("J_modal");
      xValueManagerDist<double>::map_const_iterator itm = double_manager.begin();
      xValueManagerDist<double>::map_const_iterator itmEnd = double_manager.end();

      std::map<string, vector<double>> res;

      for (; itm != itmEnd; ++itm)
      {
         assert(test_1d(itm->first.getPhys()));
         string mode_name = xKeyInfo::getGeomName(itm->first.getGeom());

         xcValueJsAndSifs::sifs(1);
         double K1 = itm->second->getVal();
         xcValueJsAndSifs::sifs(2);
         double K2 = itm->second->getVal();
         xcValueJsAndSifs::sifs(3);
         double K3 = itm->second->getVal();
         xcValueJsAndSifs::js(1);
         double J = itm->second->getVal();

         std::vector<double> v;
         v.push_back(K1);
         v.push_back(K2);
         v.push_back(K3);
         v.push_back(J);
         res[mode_name] = v;
      }

      std::map<string, vector<double>>::const_iterator its = res.begin();
      std::map<string, vector<double>>::const_iterator itse = res.end();
      for (; its != itse; ++its)
      {
         string s = its->first;
         vector<double> v = its->second;
         out.width(12);
         out << setprecision(3) << s;
         for (size_t i = 0; i < v.size(); ++i)
         {
            out.width(12);
            out << v[i];
         }
         out << endl;
      }
      out.close();
   }
}

// void CreateSubsetforJint(xMesh & mesh, const xMesh &front_mesh, std::function<double (mVertex*)> front_distance,
// 			 double rho_cylinder, int nb_layers_cylinder,
// 			 int nb_layers_core, const std::string &subset_for_q, const std::string &subset_for_int ){

//   mesh.allocateSubsetEntities("cylinder_geo_around_front");
//   int dim = mesh.dim();
//   if (mesh.size(0) == 0) std::cout << "warning : the mesh is empty" << std::endl;
//   for(xIter it=mesh.begin(0); it!=mesh.end(0); ++it)
//     {
//       mVertex *v  = (mVertex*) *it;
//       double r = front_distance(v);
//       if (r <= rho_cylinder)
// 	{
// 	  for (int i=0; i < v->size(dim); ++i)
// 	    {
// 	      mEntity* e = v->get(dim,i);
// 	      mesh.add_sub(e, "cylinder_geo_around_crack_front");
// 	      for (int ii=0; ii < dim; ++ii)
// 		{
// 		  for (int j=0; j < e->size(ii); ++j)
// 		    {
// 		      mEntity* vv = e->get(ii,j);
// 		      mesh.add_sub(vv, "cylinder_geo_around_front");
// 		    }
// 		}
// 	    }
// 	}
//     }
//     if (mesh.size_sub(dim, "cylinder_geo_around_front") == 0)
//     {
//       cout << " You must increase the size of rho for the domain integral because the domain is empty " << endl;
//       cout << " Using topologic definition to build a domain " << endl;
//     }

//     xcElementsAlongFrontCreator seed_creator(front_mesh, front_distance);
//     mesh.createSubsetEntities("cylinder_core_seed", seed_creator);

//     xAddLayerCreator add_layer_creator(nb_layers_cylinder+nb_layers_core, "cylinder_core_seed", 0);
//     mesh.createSubsetEntities("cylinder_topo_around_front", add_layer_creator);

//      // domain for integral if the union of the geo and topo cylinder
//   // around the front MINUS the core

//   // core definition
//   xAddLayerCreator add_layer_creator_for_core(nb_layers_core, "cylinder_core_seed", 0);
//   mesh.createSubsetEntities("cylinder_core", add_layer_creator_for_core);

//   //union
//   xUnionCreator union_creator("cylinder_geo_around_front", "cylinder_topo_around_front");
//   mesh.createSubsetEntities(subset_for_q, union_creator);

//   //cylinder minus core : domain for int
//   xFirstMinusSecondCreator minus_creator(subset_for_q, "cylinder_core");
//   mesh.createSubsetEntities(subset_for_int, minus_creator);

//   mesh.modifyState_sub(dim,0,true, subset_for_int);
//   mesh.modifyState_sub(dim,dim-1,true, subset_for_int);
//   mesh.modifyState_sub(dim-1,dim,true, subset_for_int);

//   // cleaning up the temp subset
//   mesh.removeSubsetEntities("cylinder_geo_around_front");
//   mesh.removeSubsetEntities("cylinder_core_seed");
//   mesh.removeSubsetEntities("cylinder_topo_around_front");
//   mesh.removeSubsetEntities("cylinder_core");
// }

// void CreateSubsetforJint(xMesh & mesh, const lCrack &crack, const mVertex &vseed,  double rho_cylinder, int nb_layers_cylinder,
// 			 int nb_layers_core, const std::string &subset_for_q, const std::string &subset_for_int ){
//  std::cout << "createSubSetforJint V2 " << std::endl;
//  mesh.allocateSubsetEntities("cylinder_geo_around_crack_front");
//  int dim = crack.dim();
//  if (dim != 2) throw;
//  xtensor::xPoint pseed = vseed.point();
//  //std::cout << pseed(0) << " " << pseed(1) << std::endl;
//  mesh.modifyState(0, 2 , true );
//  for(xIter it=mesh.begin(0); it!=mesh.end(0); ++it)
//     {
//       mEntity *v  = *it;
//       xtensor::xPoint pcourant = ((mVertex *)(v))->point();
//       double a=(pcourant(0)-pseed(0));
//       double b = (pcourant(1)-pseed(1));
//       double r = sqrt(a*a + b*b);
//       // std::cout << r << " " << rho_cylinder << std::endl;
//       if (r <= rho_cylinder)
// 	{
// 	  for (int i=0; i < v->size(dim); ++i)
// 	    {
// 	      mEntity* e = v->get(dim,i);
// 	      mesh.add_sub(e, "cylinder_geo_around_crack_front");
// 	      for (int ii=0; ii < dim; ++ii)
// 		{
// 		  for (int j=0; j < e->size(ii); ++j)
// 		    {
// 		      mEntity* vv = e->get(ii,j);
// 		      mesh.add_sub(vv, "cylinder_geo_around_crack_front");
// 		    }
// 		}
// 	    }
// 	}
//     }
//  mesh.export_sub( "cylinder_geo_around_crack_front.msh", "cylinder_geo_around_crack_front");
//     if (mesh.size_sub(dim, "cylinder_geo_around_crack_front") == 0)
//     {
//       cout << " You must increase the size of rho for the domain integral because the domain is empty " << endl;
//       cout << " Using topologic definition to build a domain " << endl;
//     }

//     xUpperCreatorRecursive seedtoelem(2);

//     mEntity *elem = seedtoelem(const_cast<mVertex * > (&vseed));
//     std::cout << "bla" << elem  << std::endl;
//     mesh.allocateSubsetEntities("cylinder_core_seed");
//     mesh.add_sub(elem,"cylinder_core_seed");

//     mesh.allocateSubsetEntities("cylinder_topo_around_crack_front");
//     xAddLayerCreator add_layer_creator(nb_layers_cylinder+nb_layers_core, "cylinder_core_seed", 0);
//     mesh.createSubsetEntities("cylinder_topo_around_crack_front", add_layer_creator);

//      // domain for integral if the union of the geo and topo cylinder
//   // around the front MINUS the core

//   // core definition
//   xAddLayerCreator add_layer_creator_for_core(nb_layers_core, "cylinder_core_seed", 0);
//   mesh.createSubsetEntities("cylinder_core", add_layer_creator_for_core);

//   //union

//   xUnionCreator union_creator("cylinder_geo_around_crack_front", "cylinder_topo_around_crack_front");
//   mesh.createSubsetEntities(subset_for_q, union_creator);

//   //cylinder minus core : domain for int
//   xFirstMinusSecondCreator minus_creator(subset_for_q, "cylinder_core");
//   mesh.createSubsetEntities(subset_for_int, minus_creator);

//   //mesh.export_sub( "cylinder_geo_around_crack_front.msh", "cylinder_geo_around_crack_front");
//   //mesh.export_sub( "cylinder_topo_around_crack_front.msh", "cylinder_topo_around_crack_front");

//   mesh.modifyState_sub(dim,0,true, subset_for_int);
//   mesh.modifyState_sub(dim,dim-1,true, subset_for_int);
//   mesh.modifyState_sub(dim-1,dim,true, subset_for_int);

//   // cleaning up the temp subset
//   mesh.removeSubsetEntities("cylinder_geo_around_crack_front");
//   mesh.removeSubsetEntities("cylinder_core_seed");
//   mesh.removeSubsetEntities("cylinder_topo_around_crack_front");
//   mesh.removeSubsetEntities("cylinder_core");

// }

// that the one

/*



  if (dim()!=2)
      printSifsModal(filename, J_modal_l, 0);
    else
      printSifsModal(filename, J_modal_l, frontseed);

    if (dim()==2){
      xEvalField<xtool::xIdentity<double> > eval(J_modal_l);
      xtensor::xPoint xyz = ((mVertex * )frontseed)->point();
      double s =0.;
      xfem::xGeomElem geom_integ(frontseed);
      xUpperCreatorRecursive integ2appro(dim());
      mEntity* e_appro = integ2appro(frontseed);
      xfem::xGeomElem geom_appro(e_appro);
      geom_integ.setUVW ({0,0,0} );
      if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
      double val_j, val_k1, val_k2, val_k3;
      xcValueJsAndSifs::js(1);
      eval(&geom_appro, &geom_integ, val_j);
      xcValueJsAndSifs::sifs(1);
      eval(&geom_appro, &geom_integ, val_k1);
      xcValueJsAndSifs::sifs(2);
      eval(&geom_appro, &geom_integ, val_k2);
      //xcValueJsAndSifs::sifs(3);
      //eval(&geom_appro, &geom_integ, val_k3);
      double theta = xcrack::AngleWithMaxHoopStress(val_k1, val_k2);
      xtensor::xPoint local;
      xtensor::xVector<> eo1, eo2, eo3;
      getLocalSmoothOrthoAxis(e_appro, geom_appro.getUVW(), local, eo1, eo2, eo3);
      //v1D(frontseed) = 0.001*val_j*(eo1*cos(theta)+eo2*sin(theta));
      v1D(frontseed) = 0.03*(eo1*cos(theta)+eo2*sin(theta));
      std::cout << "res front " <<  val_j <<" " << val_k1 << " " << val_k2 << " " <<  theta << " " << eo1 << " " <<eo2  <<" " <<
std::endl;
    }
    else {
      xEvalField<xtool::xIdentity<double> > eval(J_modal_l);
      const xValueManagerDist<double>& double_manager = *J_modal_l.getValueManager();
      xIsPhys test_1d("J_modal");
      std::map<string, double > resJ;
      xValueManagerDist<double>::map_const_iterator itm = double_manager.begin();
      for ( ; itm != double_manager.end(); ++itm) {
        assert(test_1d(itm->first.getPhys()));
        string mode_name = xKeyInfo::getGeomName(itm->first.getGeom());
        xcValueJsAndSifs::js(1);
        double J = itm->second->getVal();
        resJ[mode_name] = J;
      }

      double Jmid = resJ.begin()->second;
      string mode_name = resJ.begin()->first;
      std::cout << mode_name <<" " << Jmid  <<std::endl;
      //throw;
      xIter it=  mesh_crack_front->begin(0);

      for (; it !=  mesh_crack_front->end(0); ++it){
        mVertex * frontvertex = (mVertex * ) (*it);
    xtensor::xPoint xyz = frontvertex->point();
        double s =0.;
        xfem::xGeomElem geom_integ(frontvertex);
        xUpperCreatorRecursive integ2appro(dim());
        mEntity* e_appro = integ2appro(frontvertex);
        xfem::xGeomElem geom_appro(e_appro);
    geom_integ.setUVW ({0,0,0} );
        if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
        double val_j, val_k1, val_k2, val_k3;
        xcValueJsAndSifs::js(1);
        eval(&geom_appro, &geom_integ, val_j);
        xcValueJsAndSifs::sifs(1);
        eval(&geom_appro, &geom_integ, val_k1);
        xcValueJsAndSifs::sifs(2);
        eval(&geom_appro, &geom_integ, val_k2);
        //xcValueJsAndSifs::sifs(3);
        //eval(&geom_appro, &geom_integ, val_k3);
        double theta = xcrack::AngleWithMaxHoopStress(val_k1, val_k2);
    xtensor::xPoint local;
        xtensor::xVector<> eo1, eo2, eo3;
        getLocalSmoothOrthoAxis(e_appro, geom_appro.getUVW(), local, eo1, eo2, eo3);
        //v1D(frontseed) = 0.001*val_j*(eo1*cos(theta)+eo2*sin(theta));
        v1D(frontvertex) = 0.05*val_j/Jmid*(eo1*cos(theta)+eo2*sin(theta));
        std::cout << "res front " << "Jmid " << Jmid  << "val_j "<< val_j  <<" " << "valk1" <<  val_k1 << " " << val_k2 << " " <<
theta << " " << eo1 << " " <<eo2  <<" " << std::endl;
      }
    }





    delete(modal);


  }

  return;

}*/

}  // namespace xcrack
// end namespace
// void lCrack::getJint3DModalParks(xField<>& disp_l ,
// 			         xIntegrationRule& integrator_vol,
// 			         xIntegrationRule& integrator_bnd,
// 			         const double& rho_cylinder,
// 				 int nb_layers_cylinder,
// 				 int nb_layers_core,
// 				 int nb_modes,
// 			         std::string filenamebase,
// 				 const xcEvalJintTerms& terms_evaluator
// 				 )
// {
//   //  assert(dim() == 3);
//   const bool debug = false;
//   if (debug) cout <<"Starting Parks Modal method "<<endl;

//   //check if lss have been created

//   mesh_crack_front->modifyAllState();
//   int ncrackpart =1;
//   if (dim()==2)
//     ncrackpart = mesh_crack_front->size(0);
//   xIter itseed  = mesh_crack_front->begin(0);
//   std::string filename;
//   std::cout << "Nb parts on the front " << ncrackpart << std::endl;
//   for (int icrackpart = 0; icrackpart < ncrackpart; icrackpart++){
//     if (ncrackpart == 1) filename =  filenamebase;
//     else {
//       stringstream ssfilename;
//       ssfilename <<filenamebase <<   "_part_" << icrackpart ;
//       filename = ssfilename.str();
//     }
//     //    std::cout << "FILENAME " << filename << std::endl;

//     mVertex *frontseed = (mVertex *) *itseed ;
//     itseed++;

//     //creation de l'ensemble des elements utiles pour l'integrale J et la definition de q
//     std::string subset_for_q("subset_for_q");
//     std::string subset_for_int("subset_for_int");
//     if (dim()==3)
//       CreateSubsetforJint(*mesh, *this, rho_cylinder, nb_layers_cylinder, nb_layers_core, subset_for_q,  subset_for_int );
//     else
//       CreateSubsetforJint(*mesh, *this, *frontseed, rho_cylinder, nb_layers_cylinder, nb_layers_core, subset_for_q,
//       subset_for_int );

//     xRegion domain_for_q(mesh, subset_for_q);
//     xRegion domain_for_integral(mesh, subset_for_int);

//     mesh->export_sub(filename + "_domain_for_integral.msh", subset_for_int);

//     xFilteredRegion<xIter, xAcceptOnBoundary> domain_for_integral_ext_bnd(domain_for_integral.begin(dim()-1),
// 									  domain_for_integral.end(dim()-1), xAcceptOnBoundary());
//     std::cout << "lss_created " << lss_created << std::endl;
//     std::cout << "checkIfClosed()" << checkIfClosed() << std::endl;

//     if (!lss_created)  setLssOnFront(-1.,1.);
//     //transformation of lss and lss_1D into angle informations.

//     if (debug){
//       xexport::xExportGmshAscii pexport;
//       ExportLevelSet(lss_1D, pexport, "lss_1D_resized");
//       ExportLevelSet(lss, pexport, "lss_resized");

//       if (checkIfClosed()){
// 	ExportLevelSet(lss_1D_cos, pexport, "lss_1D_cos");
// 	ExportLevelSet(lss_1D_sin, pexport, "lss_1D_sin");
// 	ExportLevelSet(lss_cos, pexport, "lss_cos");
// 	ExportLevelSet(lss_sin, pexport, "lss_sin");
//       }
//     }
//     //if (dim() != 3) throw;
//     if (dim()==2) nb_modes =1;
//     xValueManagerDist<double> double_manager;
//     xSpaceRegular * modal = 0;

//     xField<> J_modal_l(&double_manager);

//     if (checkIfClosed()){
//       modal = new xcSpaceModalFourier("J_modal", xSpace::SCALAR, nb_modes, lss_cos, lss_sin);
//       J_modal_l.insert( *((xcSpaceModalFourier *) modal) );
//     }
//     else{
//       modal = new xcSpaceModalLegendre("J_modal", xSpace::SCALAR, nb_modes, lss);
//       J_modal_l.insert( *((xcSpaceModalLegendre *) modal) );
//     }

//     //declaration du systeme

//     mEntity * e = *mesh->begin(dim());
//     xfem::xGeomElem geo(e);
//     xMaterial *mat = xMaterialManagerSingleton::instance().getMaterial(&geo);

//     const xTensors* properties = mat->getProperties();
//     double young   = properties->scalar("YOUNG_MODULUS");
//     double poisson = properties->scalar("POISSON_RATIO");
//     double Estar = young/(1.-poisson*poisson);

//     xcValueJs::setNbModes(1);
//     xcValueSifs::setNbModes(3);
//     xcValueSifs::setGeom(xcValueSifs::GEOM_3D);
//     if (dim()==2)  xcValueSifs::setGeom(xcValueSifs::GEOM_PLANE_STRAIN);
//     xcValueSifs::setYoungAndPoisson(young, poisson);

//     //Declaration de l'interpolation (xField<> + support)

//     xValueCreator<xcValueJsAndSifs>  creator_double;
//     //xValueCreator<xValueDouble > creator_double;

//     DeclareInterpolation(J_modal_l, creator_double, mesh_crack_front->begin(0), ++mesh_crack_front->begin(0));
//     xStateDofCreator<> snh(double_manager, "dofs");

//     DeclareState(J_modal_l, snh, mesh_crack_front->begin(0), ++mesh_crack_front->begin(0));

//     /*  if (debug){
//       int nddl = double_manager.size("dofs");
//       for (int kk=0; kk< nddl; ++kk){
// 	std::cout << "exporting mode " << kk << std::endl;
// 	xlinalg::xCSRVector sol(nddl);
// 	sol[kk] = 1.;
// 	Visit(xWriteSolutionVisitor(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));
// 	xEvalField<xtool::xIdentity <double> > eval_modal( J_modal_l) ;
// 	xEvalGradField<xtool::xIdentity <xtensor::xVector<> > > eval_modal_grad( J_modal_l) ;

// 	xexport::xExportGmshAscii pexport;
// 	//	xRegion all(data.mesh);
// 	//	xIntegrationRulePartition   integrator_exp(1);
// 	xIntegrationRuleBasic   integrator_exp(1);
// 	std::stringstream filenamemo;
// 	filenamemo << "modal_"<< kk;
// 	Export(eval_modal, pexport, filenamemo.str(), integrator_exp, domain_for_integral.begin(3),  domain_for_integral.end(3));
// 	filenamemo << "_grad" ;
// 	Export(eval_modal_grad, pexport, filenamemo.str(), integrator_exp, domain_for_integral.begin(3),
// domain_for_integral.end(3));
//       }
//       xlinalg::xCSRVector sol(nddl);
//       Visit(xWriteSolutionVisitor(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));
//     }
//     */

//     if (debug)
//       double_manager.PrintForDebug("dcl_J_modal.dbg");

//     //const xLevelSet& lsn = *getFieldn();
//     //const xLevelSet& lst = *getFieldt();

//     //fonction fr pour Parks
//     xLevelSet fr(domain_for_q);
//     xcSetRadialFunction set_fr;
//     fr.accept(set_fr);
//     xEvalLevelSet<xtool::xIdentity<double> > eval_fr(fr);
//     xEvalGradLevelSet<xtool::xIdentity<xtensor::xVector<> > > eval_grad_fr(fr);

//     xexport::xExportGmshAscii pexport;

//     //q_vector est la direction de grad_lst partout
//     //sauf sur les bords ou on ne garde que la partie
//     //tangente à la surface
//     xEvalGradSmoothLevelSet<xtool::xIdentity<xtensor::xVector<> > > eval_q_dir(*getFieldt());

//     xEvalHessianLevelSet <xtool::xIdentity<xtensor::xTensor2<> > > eval_grad_q_dir(*getFieldt());
//     xEvalBinary<xtool::xMult<xtensor::xTensor2<>, double, xtensor::xTensor2<> > > eval_der1(eval_grad_q_dir, eval_fr);
//     xEvalBinary<xTensorProduct> eval_der2(eval_q_dir, eval_grad_fr);
//     xEvalBinary<xtool::xAdd<xtensor::xTensor2<>, xtensor::xTensor2<>, xtensor::xTensor2<> > > eval_grad_q_global(eval_der1,
//     eval_der2);

//     xEvalBinary<xtool::xMult<xtensor::xVector<>, double, xtensor::xVector<> > > eval_q_global(eval_q_dir, eval_fr);
//     if (debug){
//       Export(eval_q_dir, pexport,"q_dir", integrator_vol, domain_for_q.begin(), domain_for_q.end());
//       Export(eval_q_global, pexport,"q_global", integrator_vol, domain_for_integral.begin(), domain_for_integral.end());
//       Export(eval_q_global, pexport,"q_global_bnd_ext", integrator_bnd,
// 	     domain_for_integral_ext_bnd.begin(),  domain_for_integral_ext_bnd.end(), xUpperAdjacency());
//       // Export(eval_q_global, pexport,"q_global_front", integrator_vol, mesh_crack_front->begin(dim()-2),
//       mesh_crack_front->end(dim()-2), xUpperCreatorRecursive(dim()));
//       //  Export(eval_grad_q_global, pexport,"grad_q_global", integrator_vol, domain_for_integral.begin(),
//       domain_for_integral.end());
//     }

//     // Allocate the linear system
//     xlinalg::xCSRVector RHS(double_manager.size("dofs"));
//     xlinalg::xCSRVector sol(double_manager.size("dofs"));
//     xlinalg::xCSRMatrix   M(double_manager.size("dofs"));

//     //Build the Solver
//     xlinalg::xLinearSystemSolverSuperLU<> solver;

//     xAssemblerBasic<> assembler_mat(M);

//     xIntegrationRuleBasic integration_rule_MAT(12); //check element size with total lenght of the crack and mode number
//     xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<double> >,  xValOperator<xtool::xIdentity<double> > > L2_bilinear;

//     Assemble(L2_bilinear, assembler_mat, integration_rule_MAT, J_modal_l, J_modal_l,
// 	     mesh_crack_front->begin(dim()-2), mesh_crack_front->end(dim()-2), xUpperCreatorRecursive(dim()),
// xUpperCreatorRecursive(dim()));

//     if (dim()==2) {
//       M.SoftZeroMatrix();
//       M.AddMatrix(1, 1, 1.);
//       // M.AddMatrix(0, 0, 1.);
//     }

//     if (debug)
//       M.OutputMatrixOctaveFormat("mass_matrix.m");

//     xAssemblerBasic<> assembler_rhs(RHS);

//     xEvalGradField<xtool::xIdentity<xtensor::xTensor2<> > > eval_grad_disp(disp_l);
//     //  xUniformMaterialSensitivity<xtensor::xTensor4<> > hooke("strain");
//     xEvalField<xtool::xIdentity<double> > eval_sol(J_modal_l);
//     xEvalBinary<xtool::xMult<xtensor::xVector<>, double, xtensor::xVector<> > >eval_sol_vec(eval_q_global, eval_sol);

//     // J computation
//     cout << " solving for J energy release "  << endl;
//     xcValueJsAndSifs::js(1);
//     const xEval<xtensor::xTensor2<> > & eval_eshelby =  terms_evaluator.get_eshelby_evaluator();
//     const xEval<xtensor::xVector<> >  & eval_source_term = terms_evaluator.get_source_evaluator();
//     xEvalNormal eval_normal;
//     xEvalBinary<xtool::xMult<xtensor::xTensor2<>, xtensor::xVector<>, xtensor::xVector<> > > eval_eshelby_n(eval_eshelby,
//     eval_normal); xEvalBinary<xtool::xMult<xtensor::xVector<>, xtensor::xVector<>, double> > eval_q_dir_eshelby_n(eval_q_dir,
//     eval_eshelby_n);

//     xFilteredRegion<xIter, xAcceptOnBoundary> domain_for_q_bnd(domain_for_q.begin(dim()-1),
// 							       domain_for_q.end(dim()-1), xAcceptOnBoundary());
//     if(debug){
//       Export(eval_q_dir_eshelby_n, pexport,"eval_q_dir_eshelby_n_ext_bnd",
// 	     integrator_bnd, domain_for_q_bnd.begin(),
// 	     domain_for_q_bnd.end(), xUpperAdjacency());
//     }

//     //volumic term
//     xcFormLinearEnergyRelease form_linear_J           (eval_q_global, eval_grad_q_global, eval_eshelby, eval_source_term);

//     Assemble(form_linear_J, assembler_rhs, integrator_vol, J_modal_l, domain_for_integral.begin(), domain_for_integral.end() );

//     //boundary term
//     xcFormLinearEnergyReleaseBoundary form_linear_interaction_boundary(eval_q_global, eval_eshelby);

//     Assemble(form_linear_interaction_boundary, assembler_rhs, integrator_bnd, J_modal_l,
//   	     domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end(),  xUpperAdjacency());
//     solver.connectMatrix(M);
//     solver.solve(sol, RHS);

//     Visit(xWriteSolutionVisitor<>(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));

//     if (debug)
//       double_manager.PrintForDebug("sol_J.dbg");

//     // Export(eval_sol_vec, pexport,"J_vec", integrator_vol, mesh_crack_front->begin(dim()-2), mesh_crack_front->end(dim()-2),
//     xUpperCreatorRecursive(dim()));

//     // Export(eval_sol,     pexport,"J_sca", integrator_vol, mesh_crack_front->begin(dim()-2), mesh_crack_front->end(dim()-2),
//     xUpperCreatorRecursive(dim()));

//     //Sifs computation
//     // mode refer to the opening mode.
//     for (int mode = 1; mode <= 3; ++mode)
//       {
// 	cout << " solving sifs for mode " << mode << endl;
// 	terms_evaluator.setmode(mode);
// 	xcValueJsAndSifs::sifs(mode);
// 	const xEval<xtensor::xTensor2<> > &  eval_eshelby_interaction = terms_evaluator.get_interaction_evaluator();
// 	const xEval<xtensor::xVector<> > & eval_source_interaction = terms_evaluator.get_source_interaction_evaluator();

// 	RHS.ZeroArray();
//         //volumic term
// 	xcFormLinearEnergyRelease form_linear_interaction(eval_q_global, eval_grad_q_global, eval_eshelby_interaction,
// eval_source_interaction ); 	Assemble(form_linear_interaction, assembler_rhs, integrator_vol, J_modal_l,
// 		 domain_for_integral.begin(), domain_for_integral.end());

// 	//boundary term
// 	xcFormLinearEnergyReleaseBoundary form_linear_interaction_boundary(eval_q_global, eval_eshelby_interaction);
// 	Assemble(form_linear_interaction_boundary, assembler_rhs, integrator_bnd, J_modal_l,
// 		 domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end(),  xUpperAdjacency());

// 	sol.ZeroArray();

// 	//WARNING THis is for K1 and K2 only ...
// 	//For K3 it should be  K3 = 2.*mu * I3 /(2.0);

// 	solver.connectMatrix(M);
// 	solver.solve(sol, RHS);
// 	transform(sol.begin(), sol.end(), sol.begin(), bind2nd(multiplies<double>(), 0.5*Estar) );

// 	Visit(xWriteSolutionVisitor<>(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));

//         std::stringstream fic_sca;
//         fic_sca << "sif_parks_mode_"  << mode << "_vec";

// 	//Export(eval_sol_vec, pexport, fic_sca.str(), integrator_vol,
// 	//       mesh_crack_front->begin(dim()-2), mesh_crack_front->end(dim()-2), xUpperCreatorRecursive(dim()) );

//         std::stringstream fic_vec;
//         fic_sca << "sif_parks_mode_"  << mode << "_vec";

// 	//Export(eval_sol, pexport, fic_vec.str(), integrator_vol,
// 	//       mesh_crack_front->begin(dim()-2), mesh_crack_front->end(dim()-2), xUpperCreatorRecursive(dim()) );

//       }

//     //Err export
//     // xcValueJsAndSifs::err();
//     //    Export(eval_sol_vec, pexport,"K_err_vec", integrator_vol, mesh_crack_front->begin(dim()-2),
//     mesh_crack_front->end(dim()-2), xUpperCreatorRecursive(dim())); if (dim()!=2)
//       printSifsModal(filename, J_modal_l, 0);
//     else
//       printSifsModal(filename, J_modal_l, frontseed);

//     if (dim()==2){
//       xEvalField<xtool::xIdentity<double> > eval(J_modal_l);
//       xtensor::xPoint xyz = ((mVertex * )frontseed)->point();
//       double s =0.;
//       xfem::xGeomElem geom_integ(frontseed);
//       xUpperCreatorRecursive integ2appro(dim());
//       mEntity* e_appro = integ2appro(frontseed);
//       xfem::xGeomElem geom_appro(e_appro);
//       geom_integ.setUVW ({0,0,0} );
//       if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
//       double val_j, val_k1, val_k2; //, val_k3;
//       xcValueJsAndSifs::js(1);
//       eval(&geom_appro, &geom_integ, val_j);
//       xcValueJsAndSifs::sifs(1);
//       eval(&geom_appro, &geom_integ, val_k1);
//       xcValueJsAndSifs::sifs(2);
//       eval(&geom_appro, &geom_integ, val_k2);
//       //xcValueJsAndSifs::sifs(3);
//       //eval(&geom_appro, &geom_integ, val_k3);
//       double theta = xcrack::AngleWithMaxHoopStress(val_k1, val_k2);
//       xtensor::xPoint local;
//       xtensor::xVector<> eo1, eo2, eo3;
//       getLocalSmoothOrthoAxis(e_appro, geom_appro.getUVW(), local, eo1, eo2, eo3);
//       //v1D(frontseed) = 0.001*val_j*(eo1*cos(theta)+eo2*sin(theta));
//       v1D(frontseed) = 0.03*(eo1*cos(theta)+eo2*sin(theta));
//       std::cout << "res front " <<  val_j <<" " << val_k1 << " " << val_k2 << " " <<  theta << " " << eo1 << " " <<eo2  <<" "
//       << std::endl;
//     }
//     else {
//       xEvalField<xtool::xIdentity<double> > eval(J_modal_l);
//       const xValueManagerDist<double>& double_manager = *J_modal_l.getValueManager();
//       xIsPhys test_1d("J_modal");
//       std::map<string, double > resJ;
//       xValueManagerDist<double>::map_const_iterator itm = double_manager.begin();
//       for ( ; itm != double_manager.end(); ++itm) {
// 	assert(test_1d(itm->first.getPhys()));
// 	string mode_name = xKeyInfo::getGeomName(itm->first.getGeom());
// 	xcValueJsAndSifs::js(1);
// 	double J = itm->second->getVal();
// 	resJ[mode_name] = J;
//       }

//       double Jmid = resJ.begin()->second;
//       string mode_name = resJ.begin()->first;
//       std::cout << mode_name <<" " << Jmid  <<std::endl;
//       //throw;
//       xIter it=  mesh_crack_front->begin(0);

//       for (; it !=  mesh_crack_front->end(0); ++it){
// 	mVertex * frontvertex = (mVertex * ) (*it);
// 	xtensor::xPoint xyz = frontvertex->point();
// 	//double s =0.;
// 	xfem::xGeomElem geom_integ(frontvertex);
// 	xUpperCreatorRecursive integ2appro(dim());
// 	mEntity* e_appro = integ2appro(frontvertex);
// 	xfem::xGeomElem geom_appro(e_appro);
// 	geom_integ.setUVW ({0,0,0} );
// 	if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
// 	double val_j, val_k1, val_k2; // val_k3;
// 	xcValueJsAndSifs::js(1);
// 	eval(&geom_appro, &geom_integ, val_j);
// 	xcValueJsAndSifs::sifs(1);
// 	eval(&geom_appro, &geom_integ, val_k1);
// 	xcValueJsAndSifs::sifs(2);
// 	eval(&geom_appro, &geom_integ, val_k2);
// 	//xcValueJsAndSifs::sifs(3);
// 	//eval(&geom_appro, &geom_integ, val_k3);
// 	double theta = xcrack::AngleWithMaxHoopStress(val_k1, val_k2);
// 	xtensor::xPoint local;
// 	xtensor::xVector<> eo1, eo2, eo3;
// 	getLocalSmoothOrthoAxis(e_appro, geom_appro.getUVW(), local, eo1, eo2, eo3);
// 	v1D(frontvertex) = 0.05*val_j/Jmid*(eo1*cos(theta)+eo2*sin(theta));
// 	std::cout << "res front " << "Jmid " << Jmid  << "val_j "<< val_j  <<" " << "valk1" <<  val_k1 << " " << val_k2 << " " <<
// theta << " " << eo1 << " " <<eo2  <<" " << std::endl;
//       }
//     }

//     mesh->removeSubsetEntities(subset_for_q);
//     mesh->removeSubsetEntities(subset_for_int);

//     delete(modal);

//   }

//   return;

// }
// void lCrack::getJint3DModalParksV2(
// 			         xIntegrationRule& integrator_vol,
// 			         xIntegrationRule& integrator_bnd,
// 			         const double& rho_cylinder,
// 				 int nb_layers_cylinder,
// 				 int nb_layers_core,
// 				 int nb_modes,
// 			         std::string filenamebase,
// 				 const xcEvalJintTerms& terms_evaluator
// 				 )
// {
//   //assert(dim() == 3);
//   const bool debug = false;
//   if (debug) cout <<"Starting Parks Modal method v2 "<<endl;
//   mesh_crack_front->modifyAllState();

//   std::map<std::string, std::pair<mVertex *, mVertex* > >  frontParts = Separate1dMeshBranch(* mesh_crack_front, this->label );
//   std::map<std::string, std::pair<mVertex *, mVertex* > >::iterator frontPartsIt = frontParts.begin();
//   std::map<std::string, std::pair<mVertex *, mVertex* > >::iterator frontPartsItEnd = frontParts.end();

//   std::cout << "Jint for ########## "  << frontParts.size() << std::endl;

//   for ( ; frontPartsIt!=frontPartsItEnd; ++frontPartsIt ){
//     std::string frontName = (*frontPartsIt).first;
//     mVertex * vFrontStart = ((*frontPartsIt).second).first;
//     mVertex * vFrontEnd = ((*frontPartsIt).second).second;
//     bool isLoop = (vFrontStart==vFrontEnd);
//     std::string filename (frontName);
//     xRegion frontPart(mesh_crack_front, frontName);
//     xLevelSet lss1d, lss1dCos, lss1dSin;
//     if (isLoop && dim()==3 ){
//       lss1dCos.setSupport(frontPart);
//       lss1dSin.setSupport(frontPart);
//       Parametrize1dMeshLoop(vFrontStart, lss1dCos, lss1dSin);
//     }
//     else if (dim()==3) {
//       lss1d.setSupport(frontPart);
//       Parametrize1dMeshBranch(vFrontStart, vFrontEnd, lss1d);
//     }
//     else{
//      lss1d.setSupport(frontPart);
//      lss1d(*(frontPart.begin(0))) = 1.;
//     }
//     std::cout << "front parametrised" << std::endl;

//     //creation de l'ensemble des elements utiles pour l'integrale J et la definition de q
//     std::string subset_for_q("subset_for_q");
//     std::string subset_for_int("subset_for_int");

//     CreateSubsetforJint(*this, frontPart, rho_cylinder, nb_layers_cylinder, nb_layers_core, subset_for_q,  subset_for_int );

//     xRegion domain_for_q(mesh, subset_for_q);
//     xRegion domain_for_integral(mesh, subset_for_int);
//     xFilteredRegion<xIter, xAcceptOnBoundary> domain_for_integral_ext_bnd(domain_for_integral.begin(dim()-1),
// 									  domain_for_integral.end(dim()-1),
// 									  xAcceptOnBoundary());

//     mesh->export_sub(filename + "_domain_for_integral.msh", subset_for_int);

//     xValueManagerDist<double> double_manager;
//     xSpaceRegular * modal = 0;
//     xField<> J_modal_l(&double_manager);

//     //set lss3d
//     xLevelSet lss3d, lss3dCos, lss3dSin;
//     if (isLoop && dim()==3 ){
//       lss3dCos.setSupport( domain_for_q );
//       xVelocityInitialization lsscos_i(lss1dCos, domain_for_q);
//       lss3dCos.accept( lsscos_i, domain_for_q);

//       lss3dSin.setSupport( domain_for_q );
//       xVelocityInitialization lsssin_i(lss1dSin,  domain_for_q);
//       lss3dSin.accept( lsssin_i, domain_for_q);

//       modal = new xcSpaceModalFourier("J_modal", xSpace::SCALAR, nb_modes, lss3dCos, lss3dSin);
//       J_modal_l.insert( *((xcSpaceModalFourier *) modal) );
//     }
//     else if (dim()==3) {
//       lss3d.setSupport(  domain_for_q );
//       xVelocityInitialization lss_i(lss1d,  domain_for_q);
//       lss3d.accept( lss_i, domain_for_q);
//       modal = new xcSpaceModalLegendre("J_modal", xSpace::SCALAR, nb_modes, lss3d);
//       J_modal_l.insert( *((xcSpaceModalLegendre *) modal) );
//     }
//     else {
//       lss3d.setSupport(  domain_for_q );
//       xVelocityInitialization lss_i(lss1d,  domain_for_q);
//       lss3d.accept( lss_i, domain_for_q);
//       modal = new xcSpaceModalLegendre("J_modal", xSpace::SCALAR, 1., lss3d);
//       J_modal_l.insert( *((xcSpaceModalLegendre *) modal) );
//     }

//     //set Material properties
//     mEntity * e = *mesh->begin(dim());
//     xfem::xGeomElem geo(e);
//     xMaterial *mat = xMaterialManagerSingleton::instance().getMaterial(&geo);

//     const xTensors* properties = mat->getProperties();
//     double young   = properties->scalar("YOUNG_MODULUS");
//     double poisson = properties->scalar("POISSON_RATIO");
//     double Estar = young/(1.-poisson*poisson);
//     xcValueSifs::setGeom(xcValueSifs::GEOM_3D);
//     if (dim()==2)  xcValueSifs::setGeom(xcValueSifs::GEOM_PLANE_STRAIN);
//     xcValueSifs::setYoungAndPoisson(young, poisson); xValueCreator<xcValueJsAndSifs>  creator_double;

//     //set Interpolation
//     xcValueJs::setNbModes(1);
//     xcValueSifs::setNbModes(3);

//     DeclareInterpolation(J_modal_l, creator_double, frontPart.begin(0), ++frontPart.begin(0));
//     xStateDofCreator<> snh(double_manager, "dofs");
//     DeclareState(J_modal_l, snh, frontPart.begin(0), ++frontPart.begin(0));

//        //fonction fr for Parks
//     xLevelSet fr(domain_for_q);
//     xcSetRadialFunction set_fr;
//     fr.accept(set_fr);
//     xEvalLevelSet<xtool::xIdentity<double> > eval_fr(fr);
//     xEvalGradLevelSet<xtool::xIdentity<xtensor::xVector<> > > eval_grad_fr(fr);
//     xexport::xExportGmshAscii pexport;
//     //q_vector est la direction de grad_lst partout
//     //sauf sur les bords ou on ne garde que la partie
//     //tangente à la surface
//     xEvalGradSmoothLevelSet<xtool::xIdentity<xtensor::xVector<> > > eval_q_dir(*getFieldt());

//     xEvalHessianLevelSet <xtool::xIdentity<xtensor::xTensor2<> > > eval_grad_q_dir(*getFieldt());
//     xEvalBinary<xtool::xMult<xtensor::xTensor2<>, double, xtensor::xTensor2<> > > eval_der1(eval_grad_q_dir, eval_fr);
//     xEvalBinary<xTensorProduct> eval_der2(eval_q_dir, eval_grad_fr);
//     xEvalBinary<xtool::xAdd<xtensor::xTensor2<>, xtensor::xTensor2<>, xtensor::xTensor2<> > > eval_grad_q_global(eval_der1,
//     eval_der2);

//     xEvalBinary<xtool::xMult<xtensor::xVector<>, double, xtensor::xVector<> > > eval_q_global(eval_q_dir, eval_fr);
//     if (debug){
//       Export(eval_q_dir, pexport,"q_dir", integrator_vol, domain_for_q.begin(), domain_for_q.end());
//       Export(eval_q_global, pexport,"q_global", integrator_vol, domain_for_integral.begin(), domain_for_integral.end());
//       Export(eval_q_global, pexport,"q_global_bnd_ext", integrator_bnd,
// 	     domain_for_integral_ext_bnd.begin(),  domain_for_integral_ext_bnd.end(), xUpperAdjacency());
//     }

//     // Allocate the linear system
//     xlinalg::xCSRVector RHS(double_manager.size("dofs"));
//     xlinalg::xCSRVector sol(double_manager.size("dofs"));
//     xlinalg::xCSRMatrix   M(double_manager.size("dofs"));

//     //Build the Solver
//     xlinalg::xLinearSystemSolverSuperLU<> solver;

//     xAssemblerBasic<> assembler_mat(M);

//     xIntegrationRuleBasic integration_rule_MAT(6); //check element size with total lenght of the crack and mode number
//     xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<double> >,  xValOperator<xtool::xIdentity<double> > > L2_bilinear;
//     Assemble(L2_bilinear, assembler_mat, integration_rule_MAT, J_modal_l, J_modal_l,
// 	     mesh_crack_front->begin(dim()-2), mesh_crack_front->end(dim()-2), xUpperCreatorRecursive(dim()),
// xUpperCreatorRecursive(dim()));
//     // wARNInG PROBABLY NOT TRUE (CHANGE MESH CRACK FRONT TO SOMETHIG ELSE ...)
//     if (dim()==2) {
//       M.SoftZeroMatrix();
//       M.AddMatrix(1, 1, 1.);
//     }

//     if (debug)
//       M.OutputMatrixOctaveFormat("mass_matrix.m");
//     xAssemblerBasic<> assembler_rhs(RHS);

//     //xEvalGradField<xtool::xIdentity<xtensor::xTensor2<> > > eval_grad_disp(disp_l);
//     //  xUniformMaterialSensitivity<xtensor::xTensor4<> > hooke("strain");
//     xEvalField<xtool::xIdentity<double> > eval_sol(J_modal_l);
//     xEvalBinary<xtool::xMult<xtensor::xVector<>, double, xtensor::xVector<> > >eval_sol_vec(eval_q_global, eval_sol);

//     // J computation
//     cout << " solving for J energy release "  << endl;
//     xcValueJsAndSifs::js(1);
//     const xEval<xtensor::xTensor2<> > & eval_eshelby =  terms_evaluator.get_eshelby_evaluator();
//     const xEval<xtensor::xVector<> >  & eval_source_term = terms_evaluator.get_source_evaluator();
//     xEvalNormal eval_normal;
//     xEvalBinary<xtool::xMult<xtensor::xTensor2<>, xtensor::xVector<>, xtensor::xVector<> > > eval_eshelby_n(eval_eshelby,
//     eval_normal); xEvalBinary<xtool::xMult<xtensor::xVector<>, xtensor::xVector<>, double> > eval_q_dir_eshelby_n(eval_q_dir,
//     eval_eshelby_n);

//     xFilteredRegion<xIter, xAcceptOnBoundary> domain_for_q_bnd(domain_for_q.begin(dim()-1),
// 							       domain_for_q.end(dim()-1), xAcceptOnBoundary());
//     if(debug){
//       Export(eval_q_dir_eshelby_n, pexport,"eval_q_dir_eshelby_n_ext_bnd",
// 	     integrator_bnd, domain_for_q_bnd.begin(),
// 	     domain_for_q_bnd.end(), xUpperAdjacency());
//     }

//     //volumic term
//     xcFormLinearEnergyRelease form_linear_J           (eval_q_global, eval_grad_q_global, eval_eshelby, eval_source_term);

//     Assemble(form_linear_J, assembler_rhs, integrator_vol, J_modal_l, domain_for_integral.begin(), domain_for_integral.end() );

//     //boundary term
//     xcFormLinearEnergyReleaseBoundary form_linear_interaction_boundary(eval_q_global, eval_eshelby);

//     Assemble(form_linear_interaction_boundary, assembler_rhs, integrator_bnd, J_modal_l,
//   	     domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end(),  xUpperAdjacency());

//     solver.connectMatrix(M);
//     solver.solve(sol, RHS);

//     Visit (xWriteSolutionVisitor < > (sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));

//     if (debug)
//       double_manager.PrintForDebug("sol_J.dbg");

//         //Sifs computation
//     // mode refer to the opening mode.
//     for (int mode = 1; mode <= 3; ++mode)
//       {
// 	cout << " solving sifs for mode " << mode << endl;
// 	terms_evaluator.setmode(mode);
// 	xcValueJsAndSifs::sifs(mode);
// 	const xEval<xtensor::xTensor2<> > &  eval_eshelby_interaction = terms_evaluator.get_interaction_evaluator();
// 	const xEval<xtensor::xVector<> > & eval_source_interaction = terms_evaluator.get_source_interaction_evaluator();

// 	RHS.ZeroArray();
//         //volumic term
// 	xcFormLinearEnergyRelease form_linear_interaction(eval_q_global, eval_grad_q_global, eval_eshelby_interaction,
// eval_source_interaction ); 	Assemble(form_linear_interaction, assembler_rhs, integrator_vol, J_modal_l,
// 		 domain_for_integral.begin(), domain_for_integral.end());

// 	//boundary term
// 	xcFormLinearEnergyReleaseBoundary form_linear_interaction_boundary(eval_q_global, eval_eshelby_interaction);
// 	Assemble(form_linear_interaction_boundary, assembler_rhs, integrator_bnd, J_modal_l,
// 		 domain_for_integral_ext_bnd.begin(), domain_for_integral_ext_bnd.end(),  xUpperAdjacency());

// 	sol.ZeroArray();
// 	solver.connectMatrix(M);
// 	solver.solve(sol, RHS);
// 	//WARNING THis is for K1 and K2 only ...
// 	//For K3 it should be  K3 = 2.*mu * I3 /(2.0);

// 	transform(sol.begin(), sol.end(), sol.begin(), bind2nd(multiplies<double>(), 0.5*Estar) );

// 	Visit(xWriteSolutionVisitor<>(sol.begin()), double_manager.begin("dofs"), double_manager.end("dofs"));

//         std::stringstream fic_sca;
//         fic_sca << "sif_parks_mode_"  << mode << "_vec";

//         std::stringstream fic_vec;
//         fic_sca << "sif_parks_mode_"  << mode << "_vec";
//       }

//     exportSifs(filename, J_modal_l, frontPart );

//     xEvalField<xtool::xIdentity<double> > eval(J_modal_l);
//     //const xValueManagerDist<double>& double_manager = *J_modal_l.getValueManager();
//     xIsPhys test_1d("J_modal");
//     std::map<string, double > resJ;
//     xValueManagerDist<double>::map_const_iterator itm = double_manager.begin();
//     for ( ; itm != double_manager.end(); ++itm) {
//       assert(test_1d(itm->first.getPhys()));
//       string mode_name = xKeyInfo::getGeomName(itm->first.getGeom());
//       xcValueJsAndSifs::js(1);
//       double J = itm->second->getVal();
//       resJ[mode_name] = J;
//     }

//     double Jmid = resJ.begin()->second;
//     string mode_name = resJ.begin()->first;
//     std::cout << mode_name <<" " << Jmid  <<std::endl;
//     //throw;
//     xIter it=  frontPart.begin(0);
//     double Jmax = 0.;
//     for (; it !=  frontPart.end(0); ++it){
//       mVertex * frontvertex = (mVertex * ) (*it);
//       xtensor::xPoint xyz = frontvertex->point();
//       //double s =0.;
//       xfem::xGeomElem geom_integ(frontvertex);
//       xUpperCreatorRecursive integ2appro(dim());
//       mEntity* e_appro = integ2appro(frontvertex);
//       xfem::xGeomElem geom_appro(e_appro);
//       geom_integ.setUVW ({0,0,0} );
//       if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
//       double val_j, val_k1, val_k2; // val_k3;
//       xcValueJsAndSifs::js(1);
//       eval(&geom_appro, &geom_integ, val_j);
//       Jmax = val_j>Jmax? val_j:Jmax;
//     }
//     it=  frontPart.begin(0);
//     for (; it !=  frontPart.end(0); ++it){
//       mVertex * frontvertex = (mVertex * ) (*it);
//       xtensor::xPoint xyz = frontvertex->point();
//       //double s =0.;
//       xfem::xGeomElem geom_integ(frontvertex);
//       xUpperCreatorRecursive integ2appro(dim());
//       mEntity* e_appro = integ2appro(frontvertex);
//       xfem::xGeomElem geom_appro(e_appro);
//       geom_integ.setUVW (xtensor::xPoint (0,0,0) );
//       if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
//       double val_j, val_k1, val_k2; // val_k3;
//       xcValueJsAndSifs::js(1);
//       eval(&geom_appro, &geom_integ, val_j);
//       xcValueJsAndSifs::sifs(1);
//       eval(&geom_appro, &geom_integ, val_k1);
//       xcValueJsAndSifs::sifs(2);
//       eval(&geom_appro, &geom_integ, val_k2);
//       //xcValueJsAndSifs::sifs(3);
//       //eval(&geom_appro, &geom_integ, val_k3);
//       double theta = xcrack::AngleWithMaxHoopStress(val_k1, val_k2);
//       xtensor::xPoint local;
//       xtensor::xVector<> eo1, eo2, eo3;
//       getLocalSmoothOrthoAxis(e_appro, geom_appro.getUVW(), local, eo1, eo2, eo3);
//       v1D(frontvertex) = 0.08*val_j/Jmax*(eo1*cos(theta)+eo2*sin(theta));
//       std::cout << "res front " << "Jmid " << Jmid << "Jmax" << Jmax << "val_j "<< val_j  <<" " << "valk1" <<  val_k1 << " " <<
//       val_k2 << " " <<  theta << " " << eo1 << " " <<eo2  <<" " << std::endl;
//     }
//     mesh->removeSubsetEntities(subset_for_q);
//     mesh->removeSubsetEntities(subset_for_int);
//     delete modal;
//   }
//   return;
// }
