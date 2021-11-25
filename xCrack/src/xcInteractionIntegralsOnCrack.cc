/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// Trellis
#include "ParUtil.h"

// xlinalg
#include "xCSRMatrix.h"
#include "xCSRVector.h"
#include "xLinearSystemSolverSuperLU.h"

// xfem
#include "xAlgorithm.h"
#include "xAssembler.h"
#include "xMaterialSensitivity.h"
#include "xMesh.h"
#include "xRegion.h"
#include "xValue.h"
#include "xValueCreators.h"

// xexport
#include "xExportGmsh.h"

// xcrack
#include "JintModalParks.h"
#include "JintParks.h"
#include "lCrack.h"
#include "xcFormLinearEnergyRelease.h"
#include "xcFront.h"
#include "xcInteractionIntegralsOnCrack.h"
#include "xcValueJs.h"

using std::cout;
using std::endl;

namespace xcrack
{
xcInteractionIntegralsOnCrack::xcInteractionIntegralsOnCrack(const xcCrackBase &_crack, const xParseData &_parameters,
                                                             const std::vector<xcRhsForInteractionsIntegrals *> &_rhss)
    : crack(_crack),
      front(_crack.getMeshCrackFront(), crack.label, _crack.getMesh(), bind1st(mem_fun(&xcCrackBase::getFrontDistance), &_crack),
            _parameters),
      verbose(false),
      parameters(_parameters),
      rhss(_rhss)

{
   dim_mesh = crack.getMesh().dim();
   mEntity *e = *(crack.getMesh().begin(dim_mesh));
   xfem::xGeomElem geo(e);
   xMaterial *mat = xMaterialManagerSingleton::instance().getMaterial(&geo);
   const xTensors *properties = mat->getProperties();
   young = properties->scalar("YOUNG_MODULUS");
   poisson = properties->scalar("POISSON_RATIO");
   Estar = young / (1. - poisson * poisson);
   //  xcValueSifs::setGeom(xcValueSifs::GEOM_3D);
   // if (dim_mesh==2)  xcValueSifs::setGeom(xcValueSifs::GEOM_PLANE_STRAIN);
   // xcValueSifs::setYoungAndPoisson(young, poisson);
};

void xcInteractionIntegralsOnCrack::computeSIFS(const xField<> &disp_l, const xIntegrationRule &integrator_vol,
                                                const xIntegrationRule &integrator_bnd)
{
   xValueCreator<xValueDouble> creator;
   xcApproxFunctionNormedGradLevelSet grad_ls(*crack.getFieldt());
   // const xMesh & mesh = crack.getMesh();
   xexport::xExportGmshAscii pexport;
   xIntegrationRuleBasic intr;

   for (xcRhsForInteractionsIntegrals *rhs : rhss)
   {
      std::string fieldname = rhs->name;
      xValKeyExtend key_modifier(fieldname);
      xSpaceVectorXFEM J_modal_l(front.getSpace(), grad_ls, key_modifier);
      xField<> *fi = new xField<>(&double_manager, J_modal_l);
      std::pair<std::string, xfem::xField<> *> stringfield(fieldname, fi);
      fields_I.insert(stringfield);
      xField<> &field = *(fields_I[fieldname]);
      DeclareInterpolation(field, creator, front.getFrontMesh().begin(0), front.getFrontMesh().end(0));
      xStateDofCreator<> snh(double_manager, fieldname);
      DeclareState(field, snh, front.getFrontMesh().begin(0), front.getFrontMesh().end(0));
      //std::cout << "space size for " << fieldname << double_manager.size(fieldname) << " " << std::endl;
   }
   xField<> &field = *((*fields_I.begin()).second);
   // WARNING -> the keys are always there ... won't work if the crack is divided in 2 parts ! //get ol one frome svn ...

   const int numdofs = double_manager.size(rhss[0]->name);
   //double_manager.PrintForDebug("fieldforsifs");

   // Bilinear Form for Mass Matrix
   xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<xtensor::xVector<>>>, xValOperator<xtool::xIdentity<xtensor::xVector<>>>>
       L2_bilinear;

   // Set Up the System
   xlinalg::xCSRVector RHS(numdofs);
   xlinalg::xCSRVector sol(numdofs);
   xlinalg::xCSRMatrix M(numdofs);
   xlinalg::xLinearSystemSolverSuperLU<> solver;
   xAssemblerBasic<> assembler_mat(M);
   xIntegrationRuleBasic integration_rule_MAT(6);
   xAcceptLocalPart<xMesh::datamanager_t> filter(xMesh::get_const_was_created_by());
   int dim_front = front.getFrontMesh().dim_global();
   xFilteredRegion<xIter, decltype(filter)> localpartfront(front.getFrontMesh().begin(dim_front),
                                                           front.getFrontMesh().end(dim_front), filter);

   // xexport::xExportGmshAscii pexport;
   xIntegrationRuleBasic irule;

   // int myrank = ParUtil::Instance()->rank();
   // int mpisize = ParUtil::Instance()->size() ;

   // Assemble Mass Matrix
   M.OutputMatrixMatlabFormat(cout);

   Assemble(L2_bilinear, assembler_mat, integration_rule_MAT, field, field, localpartfront.begin(), localpartfront.end(),
            xUpperCreatorRecursive(dim_mesh), xUpperCreatorRecursive(dim_mesh));

#ifdef PARALLEL
   int nnz = M.ReturnNumberOfNonZero();
   double *start = M.ReturnValPointer();
   std::vector<double> sendbuf(nnz);
   for (int i = 0; i < nnz; ++i)
   {
      sendbuf[i] = start[i];
      start[i] = 0.;
   }
   MPI_Allreduce(&sendbuf[0], &start[0], nnz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD /* ierr*/);
#endif

   M.OutputMatrixMatlabFormat(cout);

   xAssemblerBasic<> assembler_rhs(RHS);

   // For Each Front part, assemble RHS for J
   for (xcFrontPart *pfront_part : front)
   {
      xcFrontPart &front_part = *pfront_part;
      std::cout << rhss.size() << std::endl;
      for (xcRhsForInteractionsIntegrals *rhs : rhss)
      {
         std::string fieldname = rhs->name;
         cout << "Assembling Rhs for " << front_part.getFrontName() << " " << fieldname << endl;
         assembleRhsforFrontPart(*fields_I[fieldname], integrator_vol, integrator_bnd, front_part, assembler_rhs, rhs->P, rhs->s);

#ifdef PARALLEL
         std::vector<double> sendbuf(numdofs);
         for (int i = 0; i < numdofs; ++i)
         {
            sendbuf[i] = RHS[i];
            RHS[i] = 0.;
         }
         MPI_Allreduce(&sendbuf[0], &RHS[0], numdofs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD /* ierr*/);
#endif
         solver.connectMatrix(M);
         solver.solve(RHS, sol);
         cout << fieldname << " rhs " << RHS.GetVal(1) << endl;
         Visit(xWriteSolutionVisitor<>(sol.begin()), double_manager.begin(fieldname), double_manager.end(fieldname));
         RHS.ZeroArray();
         sol.ZeroArray();
      }

      exportSIFSOnFrontPart(front_part);
   }
}

xcInteractionIntegralsOnCrack::~xcInteractionIntegralsOnCrack(){
    /* iterator it = fronts.begin();
    iterator itend = fronts.end();
    xMesh &mesh_crack_front = (const_cast<xMesh & > (crack.getMeshCrackFront()));
    while (it!=itend){
      mesh_crack_front.removeSubsetEntities((*it)->getFrontName() );
      delete (*it);
      (*it) = 0;
      ++it;
    }
    fronts.clear();
    */
};

void xcInteractionIntegralsOnCrack::exportSIFS() { cout << "  xcInteractionIntegralsOnCrack::exportSIFS To Implement" << endl; }

void xcInteractionIntegralsOnCrack::postResults(std::map<xcFrontPart *, std::vector<resultsfrontpartpoint>> &results) const
{
   results.clear();
   xcFront::const_iterator itfront = front.begin_front_parts();
   xcFront::const_iterator itfrontend = front.end_front_parts();
   // For Each Front part, assemble RHS for J
   for (; itfront != itfrontend; ++itfront)
   {
      std::vector<resultsfrontpartpoint> resultsfront;
      xcFrontPart &front_part = *(*itfront);
      const int dimfront = (front_part.fronttype == xcFrontPart::Point) ? 0 : 1;
      const xRegion &frontRegion = front_part.getFrontRegion();
      xAcceptLocalPart<xMesh::datamanager_t> filter(xMesh::get_const_was_created_by());
      xFilteredRegion<xIter, decltype(filter)> filteredfront(frontRegion.begin(dimfront), frontRegion.end(dimfront), filter);
      xUpperCreatorRecursive integ2appro(dimfront + 2);
      for (mEntity *e_integ : filteredfront)
      {
         xfem::xGeomElem geom_integ(e_integ);
         mEntity *e_appro = integ2appro(e_integ);
         xfem::xGeomElem geom_appro(e_appro);
         for (int ipt = 0; ipt < dimfront + 1; ++ipt)
         {
            resultsfrontpartpoint resfrontpoint;
            double u = (ipt == 0) ? -1. : 1.;
            geom_integ.setUVW({u, 0, 0});
            if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
            xtensor::xPoint xyz = geom_integ.getXYZ();
            double s = 0.;
            const xcFrontPartLineOpen *lineOpen = dynamic_cast<const xcFrontPartLineOpen *>(&front_part);
            if (lineOpen)
               s = lineOpen->getLss1d().operator()(e_integ->get(0, ipt));
            else
            {
               const xcFrontPartLineClose *lineClose = dynamic_cast<const xcFrontPartLineClose *>(&front_part);
               if (lineClose)
               {
                  double coss = lineClose->getLss1dCos().operator()(e_integ->get(0, ipt));
                  double sins = lineClose->getLss1dSin().operator()(e_integ->get(0, ipt));
                  s = atan2(sins, coss);
               }
            }
            resfrontpoint.s = s;
            resfrontpoint.x = xyz(0);
            resfrontpoint.y = xyz(1);
            resfrontpoint.z = xyz(2);

            std::vector<xcRhsForInteractionsIntegrals *>::const_iterator itf = rhss.begin();
            std::vector<xcRhsForInteractionsIntegrals *>::const_iterator itfe = rhss.end();
            for (; itf != itfe; ++itf)
            {
               xtensor::xVector<> val_k;
               std::string name = (*itf)->name;
               const xField<> &field = *((*(fields_I.find(name))).second);
               xEvalField<xtool::xIdentity<xtensor::xVector<>>> eval(field);
               eval(&geom_appro, &geom_integ, val_k);
               resfrontpoint.vals[name] = val_k.mag();
            }
            resultsfront.push_back(resfrontpoint);
         }
      }
      results[&front_part] = resultsfront;
   }  // end frontpart iter
}

template<class T, class U, class V>
void exportsifs(std::ostream &out, const T &parameters, const U &rhss, const V & res){
    out.setf(ios_base::scientific);
    out.setf(ios_base::right);
    out << "# "
        << "Parameters for domain integral : " << endl;
    out << "#  "
        << "nb_modes            :" << parameters.getInt("sifs_nb_modes") << endl;
    out << "#  "
        << "rho_geo             :" << parameters.getDouble("sifs_rho_geo") << endl;
    out << "#  "
        << "nb_layers_core      :" << parameters.getInt("sifs_nb_layers_core") << endl;
    out << "#  "
        << "nb_layers_cylinder  :" << parameters.getInt("sifs_nb_layers_cylinder") << endl;
    out << "#";out.setf(ios_base::right); out.width(11);    out << "s";  out.width(12);
    out << "x"; out.width(12); out << "y";  out.width(12); out << "z";
    for (auto itf : rhss)
    {
       out.width(12);
       out << itf->name;
    }
    out << std::endl;
    for (auto its : res)
    {
       double s = its.first;
       const auto & v = its.second;
       out.width(12);
       out << setprecision(3) << s;
       for (auto vi :v) { out.width(12);  out << vi; }
       out << endl;
    }
}

void xcInteractionIntegralsOnCrack::exportSIFSOnFrontPart(const xcFrontPart &front_part) const
{
   // xEvalField<xtool::xIdentity<xtensor::xVector<> > > eval(J_field);
   std::map<double, vector<double>> res;
   int dimfront = 1;
   if (front_part.fronttype == xcFrontPart::Point) dimfront = 0;

   const xRegion &frontRegion = front_part.getFrontRegion();
   xAcceptLocalPart<xMesh::datamanager_t> filter(xMesh::get_const_was_created_by());
   xFilteredRegion<xIter, decltype(filter)> filteredfront(frontRegion.begin(dimfront), frontRegion.end(dimfront), filter);
   xUpperCreatorRecursive integ2appro(dimfront + 2);
   for (mEntity *e_integ : filteredfront)
   {
      xfem::xGeomElem geom_integ(e_integ);
      mEntity *e_appro = integ2appro(e_integ);
      xfem::xGeomElem geom_appro(e_appro);
      for (int ipt = 0; ipt < 2; ++ipt)
      {
         double u = (ipt == 0) ? -1. : 1.;
         geom_integ.setUVW({u, 0, 0});
         if (geom_appro.getEntity() != geom_integ.getEntity()) geom_appro.setUVWForXYZ(geom_integ.getXYZ());
         xtensor::xPoint xyz = geom_integ.getXYZ();
         double s = 0.;
         const xcFrontPartLineOpen *lineOpen = dynamic_cast<const xcFrontPartLineOpen *>(&front_part);
         if (lineOpen)
            s = lineOpen->getLss1d().operator()(e_integ->get(0, ipt));
         else
         {
            const xcFrontPartLineClose *lineClose = dynamic_cast<const xcFrontPartLineClose *>(&front_part);
            if (lineClose)
            {
               double coss = lineClose->getLss1dCos().operator()(e_integ->get(0, ipt));
               double sins = lineClose->getLss1dSin().operator()(e_integ->get(0, ipt));
               s = atan2(sins, coss);
            }
         }

         // lss_1D.getVal(e_integ, geom_integ.getUVW());
         // cout << "J " << val_j.mag() << endl;
         std::vector<double> vals{xyz(0), xyz(1), xyz(2)};
         for (xcRhsForInteractionsIntegrals *rhs : rhss)
         {
            xtensor::xVector<> val_k;
            const xField<> &field = *((*(fields_I.find((rhs)->name))).second);
            xEvalField<xtool::xIdentity<xtensor::xVector<>>> eval(field);
            eval(&geom_appro, &geom_integ, val_k);
            vals.push_back(val_k.mag());
         }
         res[s] = vals;
      }
   }

   exportsifs(std::cout, parameters, rhss,  res);

   std::stringstream filename_s;
   int myrank = ParUtil::Instance()->rank();
   filename_s << front_part.front_part_name << myrank << "_s.txt";
   std::ofstream out(filename_s.str().c_str());

   exportsifs(out, parameters, rhss,  res);
   /*
   out.setf(ios_base::scientific);
   out.setf(ios_base::right);
   out << "# "
       << "Parameters for domain integral : " << endl;
   out << "#  "
       << "nb_modes            :" << parameters.getInt("sifs_nb_modes") << endl;
   out << "#  "
       << "rho_geo             :" << parameters.getDouble("sifs_rho_geo") << endl;
   out << "#  "
       << "nb_layers_core      :" << parameters.getInt("sifs_nb_layers_core") << endl;
   out << "#  "
       << "nb_layers_cylinder  :" << parameters.getInt("sifs_nb_layers_cylinder") << endl;

   out << "#";
   out.setf(ios_base::right);
   out.width(11);
   out << "s";
   out.width(12);
   out << "x";
   out.width(12);
   out << "y";
   out.width(12);
   out << "z";

   for (auto itf : rhss)
   {
      out.width(12);
      out << itf->name;
   }
   out << std::endl;
   for (auto its : res)
   {
      double s = its.first;
      const auto & v = its.second;
      out.width(12);
      out << setprecision(3) << s;
      for (auto vi :v) { out.width(12);  out << vi;}
      out << endl;
   }
   out.close();
   */
}

}  // namespace xcrack

/*
{
string filename_modes = filename + "_modes.txt";
std::ofstream out(filename_modes.c_str());
out.setf(ios_base::scientific);
out.setf(ios_base::right);

out.setf(ios_base::scientific);
out.setf(ios_base::right);
out.width(12); out << "mode";
out.width(12); out << "K1";
out.width(12); out << "K2";
out.width(12); out << "K3";
out.width(12); out << "J\n";



const xValueManagerDist<double>& double_manager = *j_modal.getValueManager();

xIsPhys test_1d("J_modal");
xValueManagerDist<double>::map_const_iterator itm = double_manager->begin();
xValueManagerDist<double>::map_const_iterator itmEnd = double_manager->end();

std::map<string, vector<double> > res;

for ( ; itm != itmEnd; ++itm) {
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

std::map<string, vector<double> >::const_iterator its   = res.begin();
std::map<string, vector<double > >::const_iterator itse = res.end();
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
*/
