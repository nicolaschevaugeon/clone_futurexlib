/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#include <iostream>

// AOMD Include
//#include "mVector.h"
//#include "mTensor2.h"
//#include "mTensor4.h"
#include "mEntity.h"
#include "mVertex.h"

// XFEM include
#include "xAlgorithm.h"
#include "xCommandOnGeomElem.h"
#include "xElement.h"
#include "xField.h"
#include "xGeomElem.h"
#include "xMaterial.h"
#include "xMaterialManager.h"
#include "xMesh.h"
#include "xSpace.h"
#include "xStateOfValue.h"
#include "xTensors.h"
#include "xValue.h"
#include "xValueCreators.h"
#include "xVariabManager.h"
#include "xVectorField.h"
#include "xZone.h"

// xexport
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"

// Xcrack includes
#include "CrackPostpro.h"
#include "Jint3D.h"
#include "lCrack.h"

using namespace xfem;

namespace xcrack
{
void lCrack::getJint3D(const xField<>& disp_l, xValueManagerDist<double>* DoubleManager, xIntegrationRule& integrator,
                       double width, double radius, std::string filename, int layers) const
{
   const bool debug = false;
   //      double h1 = 2. * width;
   //      double h2 = 2. * width;
   double h3 = width;
   double hr = 2 * radius;
   //      double hr2 = radius;
   int postid = 0;
   // boucel sur les segments de la fissure
   xMesh* mesh_front1D = mesh_crack_front;

   v1D.setSupport(mesh_front1D);
   char s2[10];
   sprintf(s2, "-%d", timestep);
   filename = filename + s2;
   filename = filename + ".txt";

   FILE* OUTPUT2 = fopen(filename.c_str(), "a");
   fprintf(OUTPUT2, "x y z r1 r2 r3 K1 K2 K3 J J(K1,K2,K3) K1*(J)\n");
   fclose(OUTPUT2);

   xIter it = mesh_front1D->begin(mesh_front1D->dim());
   xIter itend = mesh_front1D->end(mesh_front1D->dim());

   for (; it != itend; ++it)
   {
      mEntity* e = *it;
      //	  if (debug)
      cout << postid << endl;
      //	  if (debug)
      cout << e->getId() << endl;
      xtensor::xPoint front_loc;
      mVertex* v0;
      mVertex* v1;
      if ((mesh_front1D->dim()) == 1)
      {
         v0 = (mVertex*)e->get(0, 0);
         v1 = (mVertex*)e->get(0, 1);
         front_loc = (v0->point() + v1->point()) * 0.5;
      }
      else if ((mesh_front1D->dim()) == 0)
      {
         v0 = (mVertex*)*it;
         front_loc = v0->point();
      }

      xtensor::xPoint local;
      xtensor::xVector<> eo1, eo2, eo3;

      mEntity* e2d = xMesh::get_const_was_created_by().at(*e);
      mEntity* e3d = xMesh::get_const_was_created_by().at(*e2d);

      if (e3d->getLevel() != mesh->dim()) e3d = e3d->get(mesh->dim(), 0);

      xElement elem(e3d);
      elem.xyz2uvw(front_loc);
      getLocalSmoothOrthoAxis(e3d, elem.getUvw(), local, eo1, eo2, eo3);

      // loop over the nodes
      xSubMesh& Jdom_elements = mesh->createSubMesh("Jdom_elements");

      xRegion Jregion(mesh, "Jdom_elements");

      set<mEntity*> vol;
      list<mEntity*> bnd;

      vol.insert(e3d);  // on commence a mettre le premier element dans le volume
      Jdom_elements.add(e3d);
      bnd.push_back(e3d);  // et dans la liste des elements frontiere
      mEntity* elist;
      while (!bnd.empty())
      {
         elist = bnd.front();                                    // prend le premier element de la liste
         bnd.pop_front();                                        // le degage
                                                                 //	    cout << "a";
         for (int i = 0; i < elist->size(mesh->dim() - 1); i++)  // boucle sur les faces de l'element
         {
            //	      cout << "b" << i ;
            mEntity* tface = elist->get(mesh->dim() - 1, i);
            mEntity* elist_o;
            for (int j = 0; j < tface->size(mesh->dim()); j++)  // boucle sur les element connectes a la face
            {
               //		cout << "c" << j;
               elist_o = tface->get(mesh->dim(), j);
               if (elist != elist_o)  // recherche l'element de l'autre cote de la face
               {
                  //		  cout << "!";
                  if (vol.find(elist_o) == vol.end())  // est il deja dans le volume ?
                  {
                     //		    cout << "&";
                     for (int k = 0; k < elist_o->size(0); k++)  // si non; boucle sur tous ses noeud
                     {
                        //		      cout << "d" << k;
                        mEntity* ev = elist_o->get(0, k);
                        mVertex* v = (mVertex*)ev;
                        xtensor::xPoint p = v->point();
                        xtensor::xVector<> pos(front_loc, p);
                        double a, b, c;
                        a = (pos * eo1);
                        b = (pos * eo2);
                        c = (pos * eo3);
                        if ((sqrt(a * a + b * b) <= hr / 2.) && (fabs(c) <= h3 / 2.))  // ce noeud est il dans la boite ?
                        {
                           //			cout << "*";
                           vol.insert(elist_o);         // on met l'element en questio dans le volume
                           Jdom_elements.add(elist_o);  // idem
                           bnd.push_back(elist_o);      // et on l'ajoute a la liste des elements frontiere
                                                        //( a la fin -> parcours horizontal avant profondeur)
                           break;
                        }
                     }
                  }
                  break;
               }
            }
         }
         //	    cout << endl;
      }
      /*

                for(xIter it2 = mesh->begin(0); it2 != mesh->end(0); ++it2)
                {
                  mVertex *v = (mVertex*) *it2;
                  mEntity *ev= *it2;
              xtensor::xPoint p = v->point();
                  //test
                  mVector pos(front_loc, p);
                  double a,b,c;
                  a=(pos * eo1);
                  b=(pos * eo2);
                  c=(pos * eo3);
                  if ((sqrt(a*a+b*b) <= hr/2.) && (fabs(c) <= h3/2.))
                  {
                    for(int i=0;i<ev->size( mesh->dim());++i)
                    {
                      mEntity* e = ev->get(mesh->dim(),i);
                      mesh->add_sub(e,"Jdom_elements");
                      if (debug) cout << "elt found" << endl;
                    }
                  }
                }*/

      // on ajoute au moins les n couches d'elements autour de la fissure, meme  _ si la boite est trop petite...
      // attention  : si le front de fissure est dans un "coin"   ( |/ ) du domaine ca ne marchera de toute facon pas.

      int layer = 0;
      set<mEntity*> elems;
      set<mEntity*> bound, bound1;
      set<mEntity*>::iterator itb;
      elems.insert(e3d);
      bound.insert(e3d);

      while ((layer <= layers) && (bound.size() != 0))
      {
         for (itb = bound.begin(); itb != bound.end(); itb++)
         {
            mEntity* el = *itb;
            Jdom_elements.add(el);
            for (int i = 0; i < el->size(0); ++i)
            {
               mEntity* ev = el->get(0, i);
               for (int j = 0; j < ev->size(mesh->dim()); ++j)
               {
                  mEntity* e = ev->get(mesh->dim(), j);
                  if (!(elems.insert(e).second))
                  {
                     bound1.insert(e);
                     //		    mesh->add_sub(e,"Jdom_elements");
                  }
               }
            }
         }
         bound = bound1;
         layer++;
      }
      /*
                for(int i=0;i<e3d->size(0);++i)
                {
                  mEntity *ev= e3d->get(0,i);
                  for(int j=0;j<ev->size(mesh->dim());++j)
                    {
                      mEntity* e = ev->get(mesh->dim(),j);
                      mesh->add_sub(e,"Jdom_elements");
                      if (debug) cout << "elt found" << endl;
                    }
                }
      */
      xSubMesh& Jdom_bnd_elements = mesh->createSubMesh("Jdom_bnd_elements");
      xRegion Jregion_bnd(mesh, "Jdom_bnd_elements");
      for (xIter it3 = Jdom_elements.begin(mesh->dim()); it3 != Jdom_elements.end(mesh->dim()); ++it3)
      {
         mEntity* e = *it3;
         for (int i = 0; i < e->size(mesh->dim() - 1); ++i)
         {
            mEntity* f = e->get(mesh->dim() - 1, i);
            if ((f->size(mesh->dim()) == 1) || (Jdom_elements.find(f->get(mesh->dim(), 0)) == nullptr) ||
                (Jdom_elements.find(f->get(mesh->dim(), 1)) == nullptr))
            {
               Jdom_bnd_elements.add(f);
               if (debug) cout << "f found" << endl;
            }
         }
      }

      xSpaceLagrange lagScalar("FCTQ", xSpace::SCALAR, xSpaceLagrange::DEGREE_ONE);
      // Declaration de la fct Alpha (norme de la fct poids)
      const xField<> alpha(DoubleManager, lagScalar);

      xValueCreator<xValueDouble> creator;
      if (debug) cout << "before declare alpha" << endl;
      DeclareInterpolation(alpha, creator, Jregion.begin(), Jregion.end());

      if (debug) cout << "before volume " << endl;
      DirichletBoundaryCondition(alpha, "FCTQ", Jregion.begin(), Jregion.end(), 1.0);
      if (debug) DoubleManager->PrintForDebug("fctq-tot.dbg");

      if (debug) cout << "before bnd " << endl;
      DirichletBoundaryCondition(alpha, "FCTQ", Jregion_bnd.begin(mesh->dim() - 1), Jregion_bnd.end(mesh->dim() - 1), 0.0);
      if (debug) DoubleManager->PrintForDebug("fctq-bnd.dbg");
      xexport::xExportGmshAscii pexport;
      xEvalField<xtool::xIdentity<double>> eval_alpha(alpha);
      xIntegrationRuleBasic integ_basic;
      {
         string s = "alpha";
         char s2[24];
         sprintf(s2, "-%d-%d", timestep, postid);
         s = s + s2;
         xexport::Export(eval_alpha, pexport, s, integ_basic, Jregion.begin(), Jregion.end());
      }

      // we compute J, K1, K2, K3,

      // for debug
      xField<> J_info_l(DoubleManager, xSpaceConstant("J_INFO"));
      xValueCreator<xValueError> creator_info;
      DeclareInterpolation(J_info_l, creator_info, Jregion.begin(), Jregion.end());
      xValueError::choice("ENG");

      // we compute J, K1, K2, K3,
      xcrack::J3DCommand_c jcommand(alpha, disp_l, *this, J_info_l);
      ApplyCommandOnIntegrationRule(jcommand, integrator, Jregion.begin(), Jregion.end());
      xTensors tensor = jcommand.getTensors();

      if (debug)
      {
         // export debug
         xEvalField<xtool::xIdentity<double>> val_J_info(J_info_l);
         xIntegrationRuleBasic integrator_export(0);
         {
            string s = "lval";
            char s2[10];
            sprintf(s2, "%d", postid);
            cout << s2 << endl;
            s = s + s2;
            xexport::Export(val_J_info, pexport, s, integrator_export, Jregion.begin(), Jregion.end());
         }
      }

      // integrate alpha over the front
      // IntegrateEvalCommand_c toto;
      double denominator{0.};
      if (mesh_front1D->dim() == 1)
      {
         xIntegrateEvalCommand<xEvalField<xtool::xIdentity<double>>> integ_command(eval_alpha, denominator);
         xIntegrationRuleBasic integrator_alpha(0);
         xFilteredRegion<xIter, xcrack::xAcceptFrontInBox> front_in_box(mesh_front1D->begin(1), mesh_front1D->end(1),
                                                                        xcrack::xAcceptFrontInBox(Jregion));
         ApplyCommandOnIntegrationRule(integ_command, integrator_alpha, front_in_box.begin(), front_in_box.end(),
                                       xcrack::xFind3DAroundFront<>());
         //	      denominator = integ_command.result();
      }
      else if (mesh_front1D->dim() == 0)
      {
         xfem::xGeomElem geom_integ(e3d);
         geom_integ.setUVWForXYZ(front_loc);
         eval_alpha(&geom_integ, &geom_integ, denominator);
      }

      // devide the results by the denominator
      tensor.scalar("J") /= denominator;
      tensor.scalar("I1") /= denominator;
      tensor.scalar("I2") /= denominator;
      tensor.scalar("I3") /= denominator;
      tensor.scalar("I1aux") /= denominator;
      tensor.scalar("I2aux") /= denominator;
      tensor.scalar("I3aux") /= denominator;

      DeleteInterpolation(alpha, Jregion.begin(), Jregion.end());

      mesh->deleteSubMesh("Jdom_elements");
      mesh->deleteSubMesh("Jdom_bnd_elements");

      //    cout << "Integrale Jdom=" << tensor.scalar("J") <<endl;
      double MatInfo[2];
      MatInfo[0] = tensor.scalar("E");
      MatInfo[1] = tensor.scalar("NU");
      printf("***********************************************************************\n");
      printf("Crack is 1, Front is 1 and PostId is %d \n", postid);
      printf("Point %d on crack front: %15.7e %15.7e %15.7e\n", postid, front_loc(0), front_loc(1), front_loc(2));
      printf("MatInfo is: %15.7e %15.7e\n", MatInfo[0], MatInfo[1]);

      double Estar = MatInfo[0] / (1. - MatInfo[1] * MatInfo[1]);
      double mu = MatInfo[0] / (2. * (1. + MatInfo[1]));

      printf("J-INTEGRAL RESULT: Jdom = %22.14e\n", tensor.scalar("J"));
      double K1_modeI = sqrt(Estar * tensor.scalar("J"));
      printf("K1 from J assuming pure mode I = %22.15e\n", K1_modeI);

      double K1aux = sqrt(Estar * tensor.scalar("I1aux"));
      double K2aux = sqrt(Estar * tensor.scalar("I2aux"));
      double K3aux = sqrt(2. * mu * tensor.scalar("I3aux"));

      double K1 = Estar * tensor.scalar("I1") / (2.0);
      double K2 = Estar * tensor.scalar("I2") / (2.0);
      double K3 = 2. * mu * tensor.scalar("I3") / (2.0);
      double J_fromKs = (K1 * K1 + K2 * K2) / Estar + K3 * K3 / (2.0 * mu);

      if (debug)
      {
         printf("Interaction integral results\n");
         printf("K1  = %22.15e\n", K1);
         printf("K2  = %22.15e\n", K2);
         printf("K3  = %22.15e\n", K3);
         printf("(K1*K1 + K2*K2)/Estar + K3*K3/(2.0*mu)   = %22.15e\n", J_fromKs);
         printf("Without numerical errors, the value above should equal Jdom\n");
         printf("K1aux  = %22.15e\n", K1aux);
         printf("K2aux  = %22.15e\n", K2aux);
         printf("K3aux  = %22.15e\n", K3aux);
         printf("I11  = %22.15e\n", tensor.scalar("I1aux"));
         printf("I22  = %22.15e\n", tensor.scalar("I2aux"));
         printf("I33  = %22.15e\n", tensor.scalar("I3aux"));
      }

      FILE* OUTPUT2 = fopen(filename.c_str(), "a");
      fprintf(OUTPUT2, "%e %e %e %e %e %e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", front_loc(0), front_loc(1),
              front_loc(2), width, width, radius, K1, K2, K3, tensor.scalar("J"), J_fromKs, K1_modeI);
      Jscalar = tensor.scalar("J");
      fclose(OUTPUT2);
      postid++;

      double theta = xcrack::AngleWithMaxHoopStress(K1, K2);

      v1D(e) = (eo1 * cos(theta) + eo2 * sin(theta)) * tensor.scalar("J");
      cout << v1D(e) << endl;
   }
   if (mesh_front1D->dim() == 1)
   {
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
         cout << v1D(v) << endl;
      }
      {
         string s = "speedfront";
         char s2[10];
         sprintf(s2, "-%d", timestep);
         s = s + s2;
         v1D.exportGmsh(s);
      }
   }
}

J3DCommand_c::J3DCommand_c(const xField<>& q, const xField<>& disp, const lCrack& c, xField<>& j)
    : crack(c), eval_alpha(q), grad_alpha(q), grad_disp(disp), j_info(j)
{
   signature.register_scalar("J");
   signature.register_scalar("I1");
   signature.register_scalar("I2");
   signature.register_scalar("I3");
   signature.register_scalar("I1aux");
   signature.register_scalar("I2aux");
   signature.register_scalar("I3aux");
   signature.register_scalar("VOL");
   signature.register_scalar("E");
   signature.register_scalar("NU");
   tensors.setSignature(&signature);
   MatInfo.push_back(0.);
   MatInfo.push_back(0.);
}
// acces

void J3DCommand_c::openApproxElem(xfem::xGeomElem* g_appro) { geom_appro = g_appro; }
void J3DCommand_c::closeApproxElem(xfem::xGeomElem* g_appro)
{
   static double previous = 0.0;
   mEntity* e_appro = g_appro->getEntity();
   double c = tensors.scalar("J");
   j_info.setVal(e_appro, c - previous);
   //      total = total + c-previous;
   // cout << " total in J3 debug is " << total << endl;
   previous = c;
}

void J3DCommand_c::execute(xfem::xGeomElem* geo_integ)
{
   const bool debug = false;
   const int nb = geo_integ->GetNbIntegrationPoints();
   for (int k = 0; k < nb; k++)
   {
      geo_integ->setUVW(k);
      if (geom_appro->getEntity() != geo_integ->getEntity())
         geom_appro->setUVWForXYZ(geo_integ->getXYZ());
      else
         geom_appro->setUVW(geo_integ->getUVW());

      eval_alpha(geom_appro, geo_integ, alpha);
      grad_alpha(geom_appro, geo_integ, dalpha);

      double Jhh_l, vol_l;
      xtensor::xVector<> Ih_l;
      xtensor::xTensor2<> I_l;

      xMaterial* mat = xMaterialManagerSingleton::instance().getMaterial(geom_appro);
      if (debug) assert(mat != nullptr);
      const xTensors* properties = mat->getProperties();
      MatInfo[0] = properties->scalar("YOUNG_MODULUS");
      MatInfo[1] = properties->scalar("POISSON_RATIO");
      tensors.scalar("E") = MatInfo[0];
      tensors.scalar("NU") = MatInfo[1];
      mat->sensitivityTo("strain", DC);
      // geometrical information at the point of integration
      //
      crack.setLocalInfos(geom_appro->getEntity(), geom_appro->getUVW());
      // original calculations
      GetJint3DLevelset(MatInfo, grad_disp, alpha, dalpha, DC, geo_integ->GetWeight(), geo_integ->GetDetJac(), geom_appro, &crack,
                        Jhh_l, Ih_l, I_l, vol_l);

      tensors.scalar("J") += Jhh_l;
      tensors.scalar("I1") += Ih_l(0);
      tensors.scalar("I2") += Ih_l(1);
      tensors.scalar("I3") += Ih_l(2);
      tensors.scalar("I1aux") += I_l(0, 0);
      tensors.scalar("I2aux") += I_l(1, 1);
      tensors.scalar("I3aux") += I_l(2, 2);
      tensors.scalar("VOL") += vol_l;
   }
}

xTensors J3DCommand_c::getTensors() { return tensors; }

}  // namespace xcrack
