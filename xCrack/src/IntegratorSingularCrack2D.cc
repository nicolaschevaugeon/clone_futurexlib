/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#include <iostream>
// AOMD dep
#include "mFace.h"
#include "mVertex.h"
// xfem dep
#include "xCommandOnGeomElem.h"
#include "xElement.h"
#include "xEntityToEntity.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xMesh.h"
// xcrack dep
#include "IntegratorSingularCrack2D.h"
#include "SingularCrackMapping.h"
#include "lCrack.h"
// SolverInterface/Lapack dep
#include "xLapackInterface.h"

namespace xcrack
{
using AOMD::mEntity;
using xfem::xLevelSet;
xtensor::xVector<> IntegratorSingularCrack2D_c::orientation(mEntity* e) const
{
   AOMD::mVertex* V0 = static_cast<AOMD::mVertex*>(e->get(0, 0));
   AOMD::mVertex* V1 = static_cast<AOMD::mVertex*>(e->get(0, 1));
   AOMD::mVertex* V2 = static_cast<AOMD::mVertex*>(e->get(0, 2));
   xtensor::xVector<> vect01(V0->point(), V1->point());
   xtensor::xVector<> vect02(V0->point(), V2->point());
   return vect01 % vect02;
}

void IntegratorSingularCrack2D_c::accept(xfem::xCommandOnGeomElem& command, mEntity* e_integ) const
{
   const bool debug = false;
   if (debug)
   {
      cout << "debut integrateur singulier" << endl;
      std::cout << "On ouvre l element ";
      e_integ->print();
   }
   if (!filter_front(e_integ))
   {
      // non singular integration
      xfem::xIntegrationRulePartition integrator_partition(xfem::xMesh::getPartition, degree_gen);
      integrator_partition.accept(command, e_integ);
      return;
   }
   xfem::xPartition partition;
   xfem::xMesh::getPartition(e_integ, partition, filter);
   if (debug) std::cout << "partition size is " << partition.size() << std::endl;
   //  int dimension=e_integ->getClassification()->dim();
   const xLevelSet* lsnn = crack.getFieldn();
   const xLevelSet* lstt = crack.getFieldt();
   mEntity* elem_interp = e_integ;
   std::vector<double> valst = (lstt)->getVals(elem_interp);
   std::vector<double> valsn = (lsnn)->getVals(elem_interp);

   for (mEntity* es_integ : partition)
   {
      if (debug) std::cout << "partition found inside sub\n";
      xfem::xElement elem(e_integ);  // interpo sur element de base!
      AOMD::mVertex* Vert[3];
      AOMD::mVertex* SingVert = nullptr;
      AOMD::mVertex* s_Vert[3];
      double RAYON[3];
      // double RAYONs[3];
      bool is_node_on_crack_tip = false;
      // pour avoir la vraie orientation du geo integ
      xtensor::xVector<> orientation_es_integ = orientation(es_integ);
      //      if (debug) cout<<"dimension is : "<<e_integ->getClassification()->dim()<<endl;
      xlinalg::xDenseMatrix Aglob(4);
      xlinalg::xCSRVector b_lsn0(4);
      xlinalg::xCSRVector sol_lsn0(4);
      xlinalg::xCSRVector b_lst0(4);
      xlinalg::xCSRVector sol_lst0(4);
      xlinalg::xCSRVector v_crack_dir(3);
      for (int ii = 0; ii < 3; ii++)
      {
         Vert[ii] = (AOMD::mVertex*)e_integ->get(0, ii);
         s_Vert[ii] = (AOMD::mVertex*)es_integ->get(0, ii);
         if (debug) cout << "Vert[ii]" << Vert[ii]->point() << endl;
         elem.xyz2uvw(Vert[ii]->point());
         double val_LST = elem.getInterpoSca(valst);
         double val_LSN = elem.getInterpoSca(valsn);
         RAYON[ii] = sqrt(val_LSN * val_LSN + val_LST * val_LST);
         if (debug) cout << "RAYON[" << ii << "] " << RAYON[ii] << endl;
         //	  cout<<"RAYON["<<ii<<"] "<<RAYON[ii]<<endl;//
         Aglob(ii, 0) = (Vert[ii]->point())(0);
         Aglob(ii, 1) = (Vert[ii]->point())(1);
         Aglob(ii, 2) = (Vert[ii]->point())(2);
         Aglob(ii, 3) = 1;
         b_lsn0[ii] = val_LSN;
         b_lst0[ii] = val_LST;
         if (debug) cout << "in IntegratorSingularCrack2D.cc, val_LSN : " << val_LSN << "val_LST : " << val_LST << endl;

         /* if (dimension ==3){
         if (RAYON[ii]<1.e-14)//theoriquement=0
           {
             is_node_on_crack_tip=true;//le noeud est sur le front de fissure
             if (debug) cout<<"is_node_on_crack_tip=true"<<endl;
             SingPoint=ii;
             if (debug) cout<<"warning for this case : ii= "<<ii<<endl;
             Vert[ii]->print();
             s_Vert[ii]->print();

           }
           }*/
         int dimension = 2;
         if (dimension == 2)
         {
            if (RAYON[ii] < 1.e-14)  // theoriquement=0
            {
               is_node_on_crack_tip = true;  // le noeud est sur le front de fissure
               if (debug) cout << "is_node_on_crack_tip=true" << endl;
               //	      SingPoint=ii;
               //	      cout<<"SingPoint big = "<<SingPoint<<endl;//
               SingVert = Vert[ii];
               if (debug) cout << "SingVert = " << endl;
               if (debug) Vert[ii]->print();
            }
            if (debug) cout << "dimension is " << dimension << endl;
            if (is_node_on_crack_tip && (dimension == 2))
            {
               for (int iis = 0; iis < ii + 1; iis++)
               {
                  xtensor::xVector<> tester(SingVert->point(), s_Vert[iis]->point());
                  // cout<<"norm du tester = "<<tester.mag()<<endl;
                  if (tester.mag() < 1.e-14) SingPoint = iis;  // faire le min
               }
            }
         }
      }
      //      if (dimension == 2) //
      if (1)
      {
         Aglob(3, 0) = orientation_es_integ(0);
         Aglob(3, 1) = orientation_es_integ(1);
         Aglob(3, 2) = orientation_es_integ(2);
         Aglob(3, 3) = 0;
      }
      xlinalg::xDenseLUFactor LU(Aglob);
      solve(LU, b_lst0, sol_lst0);
      solve(LU, b_lsn0, sol_lsn0);
      for (int i = 0; i < 3; ++i)
         v_crack_dir[i] = sol_lst0[(i + 1) % 3] * sol_lsn0[(i + 2) % 3] - sol_lst0[(i + 2) % 3] * sol_lsn0[(i + 1) % 3];

      /*if (debug)
        {
          cout<<"Aglob"<<endl; Aglob.Display();
          cout<<"b_lst0"<<endl;	b_lst0.Display();
          cout<<"sol_lst0"<<endl; sol_lst0.Display();
          cout<<"b_lsn0"<<endl; b_lsn0.Display();
          cout<<"sol_lsn0"<<endl; sol_lsn0.Display();
          v_crack_dir.Display(); std::ofstream out_crack_dir("crack_dir.m");
          }*/
      xlinalg::xDenseMatrix A_get_singular_node(3);
      xlinalg::xCSRVector rhs_sing_node(3);
      xlinalg::xCSRVector sol_sing_node(3);
      for (int i = 0; i < 3; i++)
      {
         A_get_singular_node(0, i) = sol_lst0[i];
         A_get_singular_node(1, i) = sol_lsn0[i];
         A_get_singular_node(2, i) = v_crack_dir[i];
      }
      rhs_sing_node[0] = -sol_lst0[3];
      rhs_sing_node[1] = -sol_lsn0[3];
      xtensor::xPoint p0 = static_cast<AOMD::mVertex*>(es_integ->get(0, 0))->point();
      rhs_sing_node[2] = (v_crack_dir[0] * p0(0) + v_crack_dir[1] * p0(1) + v_crack_dir[2] * p0(2));
      xlinalg::xDenseLUFactor LU_sing(A_get_singular_node);
      solve(LU_sing, rhs_sing_node, sol_sing_node);
      /*
      if (debug)
        {
          cout<<"matrice pour noeud singulier"<<endl;
          A_get_singular_node.Display();
          cout<<"membre de droite"<<endl;
          rhs_sing_node.Display();
          cout<<"noeud singulier 1"<<endl;
          sol_sing_node.Display();
          std::cout << "Une partition existe\n";
          }*/

      if (is_node_on_crack_tip)
      {
         if (debug) cout << "is node on crack tip" << endl;
         Vert[SingPoint]->print();
         s_Vert[SingPoint]->print();
         SingularCrackMapping mapping(es_integ, SingPoint);
         xquadrature::xGaussIntegrator integrator(mapping.getType());
         xfem::xGeomElem geo_integ_es(es_integ, &mapping, &integrator);
         geo_integ_es.SetIntegrationPointNumberForDegree(degree);
         command.execute(&geo_integ_es);
      }
      else
      {
         // cout << " node not on ctip " << endl;
         if (debug) cout << "toto0" << endl;
         // int compteur=0;
         xfem::xMesh mesh_point;
         Trellis_Util::mPoint point_sing(sol_sing_node[0], sol_sing_node[1], sol_sing_node[2]);
         AOMD::mVertex* ve = mesh_point.getMesh().createVertex(point_sing, e_integ->getClassification());
         // debut boucle sur les edges
         for (int i = 0; i < 3; i++)
         {
            if (debug)
               cout << "points pour fabriquer la face : point 1 = " << Vert[i]->point()
                    << " point 2 = " << Vert[(i + 1) % 3]->point() << "point singulier : " << ve->point() << endl;
            if (debug) cout << "toto1" << endl;
            AOMD::mFace fac1(s_Vert[i], s_Vert[(i + 1) % 3], ve, e_integ->getClassification());
            xtensor::xVector<> vect_e1(ve->point(), s_Vert[i]->point());
            xtensor::xVector<> vect_e2(ve->point(), s_Vert[(i + 1) % 3]->point());
            xtensor::xVector<> produitvectoriel = vect_e1 % vect_e2;
            if (sqrt(produitvectoriel * produitvectoriel) > 1.e-13)
            {
               SingularCrackMapping mapping(&fac1, 2);
               xquadrature::xGaussIntegrator integrator(mapping.getType());
               xfem::xGeomElem geo_elem_superpos(&fac1, &mapping, &integrator);
               int provec_super = (produitvectoriel * orientation_es_integ > 0) ? 1 : -1;
               geo_elem_superpos.SetOrientation(provec_super);
               geo_elem_superpos.SetIntegrationPointNumberForDegree(degree);
               if (debug) cout << "degree superpos=" << degree << endl;
               if (debug)
                  cout << "nb integration points superposition(theoriquement singular)"
                       << geo_elem_superpos.GetNbIntegrationPoints() << endl;
               command.execute(&geo_elem_superpos);
            }
         }
      }
   }
   if (debug) cout << "fin integrateur singulier 2D" << endl;
}

}  // namespace xcrack
