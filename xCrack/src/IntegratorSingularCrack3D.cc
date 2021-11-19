/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#include <iostream>
#include <map>
// Trellis dep
#include "mEdge.h"
#include "mPoint.h"
#include "mTet.h"
#include "mVertex.h"
// xfem dep
#include "xCommandOnGeomElem.h"
#include "xElement.h"
#include "xEntityToEntity.h"
#include "xFemMatrix.h"
#include "xForm.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xMesh.h"
// xgeom
#include "xSimpleGeometry.h"
// xtensor
#include "xVector.h"
// xcut deps
#include "xMeshCut.h"
#include "xPhysSurf.h"
// xcrack dep
#include "IntegratorSingularCrack2D.h"
#include "IntegratorSingularCrack3D.h"
#include "SingularCrackMapping.h"
#include "lCrack.h"
// SolverInterface/Lapack dep
#include "xLapackInterface.h"

using namespace xcrack;
using namespace xfem;
using namespace xgeom;
using namespace xtensor;

xVector<> createProjectedMesh(const std::vector<xPoint> &nodes, xlinalg::xCSRVector &lsn, xlinalg::xCSRVector &lst,
                              pGEntity classification, xMesh *theMesh)
{
   xlinalg::xDenseMatrix Aglob(4);
   for (int ii = 0; ii < 4; ii++)
   {
      Aglob(ii, 0) = nodes[ii](0);
      Aglob(ii, 1) = nodes[ii](1);
      Aglob(ii, 2) = nodes[ii](2);
      Aglob(ii, 3) = 1.;
   }
   xlinalg::xCSRVector sol_lsn0(4);
   xlinalg::xCSRVector sol_lst0(4);
   xlinalg::xDenseLUFactor LU(Aglob);
   solve(LU, lst, sol_lst0);
   solve(LU, lsn, sol_lsn0);
   xlinalg::xCSRVector v_front_dir(3);
   for (int i = 0; i < 3; ++i)
      v_front_dir[i] = sol_lst0[(i + 1) % 3] * sol_lsn0[(i + 2) % 3] - sol_lst0[(i + 2) % 3] * sol_lsn0[(i + 1) % 3];
   //  v_front_dir.Cross(sol_lst0,sol_lsn0);
   double temp_double = sqrt(v_front_dir[0] * v_front_dir[0]);
   int i0 = 0;
   for (int ii = 1; ii < 3; ii++)
   {
      if (sqrt(v_front_dir[ii] * v_front_dir[ii]) > temp_double)
      {
         i0 = ii;
         temp_double = sqrt(v_front_dir[ii] * v_front_dir[ii]);
      }
   }
   xlinalg::xDenseMatrix A_get_singular_node(2);
   xlinalg::xCSRVector rhs_sing_node(2);
   xlinalg::xCSRVector sol_sing_node(2);
   // seulement 2 pour A_get singular en 3D!!! alors que 3 en 2D!
   A_get_singular_node(0, 0) = sol_lsn0[(i0 + 1) % 3];
   A_get_singular_node(0, 1) = sol_lsn0[(i0 + 2) % 3];
   rhs_sing_node[0] = -sol_lsn0[3];
   A_get_singular_node(1, 0) = sol_lst0[(i0 + 1) % 3];
   A_get_singular_node(1, 1) = sol_lst0[(i0 + 2) % 3];
   rhs_sing_node[1] = -sol_lst0[3];
   xlinalg::xDenseLUFactor LU2(A_get_singular_node);
   solve(LU2, rhs_sing_node, sol_sing_node);
   xlinalg::xCSRVector point_on_front(3);
   point_on_front[i0] = 0.;
   point_on_front[(i0 + 1) % 3] = sol_sing_node[0];
   point_on_front[(i0 + 2) % 3] = sol_sing_node[1];
   // Tous ce qui précède ne sert qu'a calculer un point sur lfront
   xPoint ve_on_front(point_on_front[0], point_on_front[1], point_on_front[2]);
   xVector<> x_front_dir(v_front_dir[0], v_front_dir[1], v_front_dir[2]);
   x_front_dir.norm();
   map<double, AOMD::mVertex *> classe;
   for (int ii = 0; ii < 4; ii++)
   {
      xVector<> P0Pii(ve_on_front, nodes[ii]);
      double absc_proj = P0Pii * x_front_dir;
      xVector<> P0Proj = x_front_dir * absc_proj;
      AOMD::mVertex *proj_vert = theMesh->getMesh().createVertex(point_on_front[0] + P0Proj[0], point_on_front[1] + P0Proj[1],
                                                                 point_on_front[2] + P0Proj[2], classification);
      classe.insert(make_pair(absc_proj, proj_vert));
   }
   map<double, AOMD::mVertex *>::iterator it = classe.begin();
   AOMD::mVertex *end = (it)->second;
   AOMD::mVertex *start;
   ++it;
   while (it != classe.end())
   {
      start = end;
      end = (it)->second;
      ++it;
      // AOMD::mEdge* e = theMesh->createEdge(start,end, classification);
      theMesh->getMesh().createEdge(start, end, classification);
   }
   return x_front_dir;
};

void IntegratorSingularCrack3D_c::accept(xCommandOnGeomElem &command3d, AOMD::mEntity *e_geom3d) const
{
   if (!filter_front(e_geom3d))
   {
      xIntegrationRulePartition integrator_partition(xMesh::getPartition, degree_gen);
      return;
   }
   xPartition partition;
   xMesh::getPartition(e_geom3d, partition, filter);
   // int dimension=e_geom3d->getClassification()->dim();
   const xLevelSet *lsnn = crack.getFieldn();
   const xLevelSet *lstt = crack.getFieldt();
   std::vector<double> valst = (lstt)->getVals(e_geom3d);
   std::vector<double> valsn = (lsnn)->getVals(e_geom3d);
   xfem::xGeomElem e_geom3d_geom(e_geom3d);
   for (AOMD::mEntity *e_sub_integ3d : partition)
   {
      xfem::xGeomElem e_geom3d_sub(e_sub_integ3d);
      xElement elem(e_geom3d);  // interpo sur element de base!
      // AOMD::mVertex* Vert[4]; -Wunused-but-set-variable says it is not used
      // AOMD::mVertex* s_Vert[4];
      std::vector<xPoint> nodes(4);
      xlinalg::xCSRVector b_lsn0(4);
      xlinalg::xCSRVector b_lst0(4);
      for (int ii = 0; ii < 4; ii++)
      {
         // Vert[ii]   = (AOMD::mVertex *)e_geom3d->get(0,ii); -Wunused-but-set-variable says it is not used
         // s_Vert[ii] = (AOMD::mVertex *)e_sub_integ3d->get(0,ii);
         nodes[ii] = ((AOMD::mVertex *)e_sub_integ3d->get(0, ii))->point();
         elem.xyz2uvw(nodes[ii]);
         b_lsn0[ii] = elem.getInterpoSca(valsn);
         b_lst0[ii] = elem.getInterpoSca(valst);
      }
      xMesh mesh_points_on_crack;
      xVector<> x_front_dir = createProjectedMesh(nodes, b_lsn0, b_lst0, e_geom3d->getClassification(), &mesh_points_on_crack);
      // if (debug)
      //	AOMD_Util::Instance()->ex_port("localcrackfront.msh",&mesh_points_on_crack);

      xMesh sub_element_3d_mesh;
      AOMD::mMesh &sub_element_3d_mmesh = sub_element_3d_mesh.getMesh();
      {
         AOMD::mVertex *tab_v[4];
         for (size_t i = 0; i < 4; i++)
         {
            Trellis_Util::mPoint p(nodes[i](0), nodes[i](1), nodes[i](2));
            tab_v[i] = sub_element_3d_mesh.getMesh().createVertex(p, e_geom3d->getClassification());
         }
         sub_element_3d_mmesh.createTetWithVertices(tab_v[0], tab_v[1], tab_v[2], tab_v[3], e_geom3d->getClassification());
         // sub_element_3d_mesh.modifyAllState();

         sub_element_3d_mmesh.modifyState(3, 2, true);
         sub_element_3d_mmesh.modifyState(3, 1, true);
         sub_element_3d_mmesh.modifyState(3, 0, true);
         sub_element_3d_mmesh.modifyState(2, 1, true);
         sub_element_3d_mmesh.modifyState(2, 0, true);
         sub_element_3d_mmesh.modifyState(1, 0, true);
         sub_element_3d_mmesh.modifyState(0, 1, true);
         sub_element_3d_mmesh.modifyState(0, 2, true);
         sub_element_3d_mmesh.modifyState(0, 3, true);
         sub_element_3d_mmesh.modifyState(1, 2, true);
         sub_element_3d_mmesh.modifyState(1, 3, true);
         sub_element_3d_mmesh.modifyState(2, 3, true);
      }

      xLevelSet lsn_sub3d, lst_sub3d;

      for (AOMD::mEntity *e : sub_element_3d_mesh.range(0))
      {
         AOMD::mVertex *Vert_test = (AOMD::mVertex *)e;
         elem.xyz2uvw(Vert_test->point());
         double val_LST2 = elem.getInterpoSca(valst);
         double val_LSN2 = elem.getInterpoSca(valsn);
         lsn_sub3d(Vert_test) = val_LSN2;
         lst_sub3d(Vert_test) = val_LST2;
      }
      // ici doit commencer la boucle sur les points de gauss 1D

      // lsn_sub3d.exportGmsh("lsn_sub3d.pos",&sub_element_3d_mesh);

      for (AOMD::mEntity *e1 : mesh_points_on_crack.range(1))
      {
         xfem::xGeomElem one_d_segment(e1);
         AOMD::mVertex *Vert00 = (AOMD::mVertex *)e1->get(0, 0);
         AOMD::mVertex *Vert01 = (AOMD::mVertex *)e1->get(0, 1);
         // si edge miniscule, i.e. probleme numerique pour la projection precedente, pas de calcul
         xVector<> seg_elem(Vert00->point(), Vert01->point());
         double norme_xvect = seg_elem.mag();
         if (norme_xvect > 1.0e-15)
         {
            int degree_on_front_direction = degree_z;
            one_d_segment.SetIntegrationPointNumberForDegree(degree_on_front_direction);
            int nb_gauss_points = one_d_segment.GetNbIntegrationPoints();
            for (int k = 0; k < nb_gauss_points; k++)
            {
               one_d_segment.setUVW(k);
               const bool special_z = false;
               xPoint integ_point;
               double w;
               if (!special_z)
                  integ_point = one_d_segment.getXYZ();
               else
               {
                  xPoint uvw = one_d_segment.getUVW();
                  w = uvw(0);
                  uvw(0) = 2. * (1. - ((1. + w) / 2.) * ((1. + w) / 2.)) - 1.;  // special Z-mapping
                  one_d_segment.setUVW(uvw);
                  integ_point = one_d_segment.getXYZ();
               }

               xPlane cutting_plane(integ_point, x_front_dir);
               xLevelSet lsortho2d(&sub_element_3d_mesh);
               lsortho2d.load(cutting_plane);
               xMesh mesh2d;
               xMesh::datamanager_t<AOMD::mEntity *> was_created_by;
               xMesh::datamanager_t<double> r_on_edge;
               xcut::createIsoZeroMeshFromLevelSet(lsortho2d, mesh2d, was_created_by, r_on_edge);
               // sub_element_3d_mesh.cutMesh(lsortho2d,&mesh2d, xtool::xIdentity<mEntity*>(),xtool::xIdentity<mEntity*>(),false);
               // // the crack surface
               xLevelSet lsn2d, lst2d;
               // mesh2d.takeTraceOn(lsn_sub3d, lsn2d);
               lsn_sub3d.takeTraceOn(mesh2d, was_created_by, r_on_edge, lsn2d);
               lst_sub3d.takeTraceOn(mesh2d, was_created_by, r_on_edge, lst2d);
               // mesh2d.takeTraceOn(lst_sub3d, lst2d);
               lCrack crack2d(&mesh2d, true);
               crack2d.setLevelSets(lsn2d, lst2d);
               IntegratorSingularCrack2D_c integrator2d(crack2d, degree, degree);  // default_filter=xAcceptAll

               xIntegrateFormCommand *f_command = dynamic_cast<xIntegrateFormCommand *>(&command3d);
               double gp_weight;
               if (!special_z)
                  gp_weight = one_d_segment.GetWeight();  // before
               else
                  gp_weight = one_d_segment.GetWeight() * (1. + w);  // special Z-mapping
               double gp_detjac = one_d_segment.GetDetJac();
               // int one_d_segment_ori=one_d_segment.GetOrientation();
               xfem::xGeomElem sub_integ3d_geom(e_sub_integ3d);
               // double gp_detjac_from_sub_integ3d_geom=sub_integ3d_geom.GetDetJac();
               xfem::xGeomElem integ3d_geom(e_geom3d);
               // double gp_detjac_from_integ3d_geom=integ3d_geom.GetDetJac();
               double coeff_accu = gp_weight * gp_detjac;
               if (f_command)
               {
                  xForm *form3d = f_command->getForm();
                  xFormBilinear<> *bf = dynamic_cast<xFormBilinear<> *>(form3d);
                  xFormLinear<> *lf = dynamic_cast<xFormLinear<> *>(form3d);
                  xFormZero<> *zf = dynamic_cast<xFormZero<> *>(form3d);
                  if (bf)
                  {
                     xFemMatrix<double> &matrix3d = bf->getMatrix();
                     matrix3d.setCoeff(coeff_accu);
                     for (AOMD::mEntity *integ2d : mesh2d.range(2)) integrator2d.accept(command3d, integ2d);
                     matrix3d.setCoeff(1.0);
                  }
                  else if (lf)
                  {
                     xFemVector<double> &vector3d = lf->getVector();
                     vector3d.setCoeff(coeff_accu);
                     for (AOMD::mEntity *integ2d : mesh2d.range(2)) integrator2d.accept(command3d, integ2d);
                     vector3d.setCoeff(1.0);
                  }
                  else if (zf)
                  {
                     xFemScalar<double> &scalar3d = zf->getScalar();
                     scalar3d.setCoeff(coeff_accu);
                     for (AOMD::mEntity *integ2d : mesh2d.range(2)) integrator2d.accept(command3d, integ2d);

                     scalar3d.setCoeff(1.0);
                  }
               }
            }
         }
      }
   }
}
