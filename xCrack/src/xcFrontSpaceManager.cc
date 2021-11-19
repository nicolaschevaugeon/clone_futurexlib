/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
// xexport
#include "xExportAlgorithm.h"
// xcrack
#include "xcFrontPart.h"
#include "xcFrontSpaceManager.h"
// xfem
#include "xParseData.h"

using namespace AOMD;
using namespace std;
using namespace xfem;

xcFrontSpaceManager::~xcFrontSpaceManager() { spaces.clear(); }

xcFrontSpaceManager::xcFrontSpaceManager(xcFrontDomainManager* _dmanager, int discretisationType, const std::string& space_name,
                                         double minimum_mode_width)
    : dmanager(_dmanager)
{
   // spaces creation
   xcFrontDomainManager::iterator it = dmanager->begin_domains_iter();
   xcFrontDomainManager::iterator iten = dmanager->end_domains_iter();
   while (it != iten)
   {
      xSpaceXFEM* space;
      const xcFrontPartBase* part = (*it)->getFrontPart();
      xValKeyExtend key_modifier("");
      switch (discretisationType)
      {
         case CONSTANT:
            (*it)->set_nb_modes(1);
            // space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR, (*(part->getFrontRegion().begin(0)))),
            // xcApproxFunctionLevelSet((*it)->getFr()), key_modifier);
            space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR, (*(part->getFrontRegion().begin(0)))),
                                   xApproxFunctionConstant(), key_modifier);
            break;
         case MODAL:  // legendre or Fourier, depending on the front type
         {
            int predefined_nb_modes = (*it)->getParameters().getInt("front_nb_modes");
            (*it)->set_nb_modes(predefined_nb_modes);
            switch (part->getFrontType())
            {
               case xcFrontPartBase::Point:
                  // space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR,
                  // (*(part->getFrontRegion().begin(0)))), xcApproxFunctionLevelSet((*it)->getFr()), key_modifier);
                  space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR, (*(part->getFrontRegion().begin(0)))),
                                         xApproxFunctionConstant(), key_modifier);
                  break;
               case xcFrontPartBase::LineOpen:
                  // space = new xSpaceXFEM(xcSpaceModalLegendre(space_name, xSpace::SCALAR, predefined_nb_modes,
                  // ((xcFrontDomainLineOpen*)(*it))->getLss3d()), xcApproxFunctionLevelSet((*it)->getFr()), key_modifier);
                  space = new xSpaceXFEM(xcSpaceModalLegendre(space_name, xSpace::SCALAR, predefined_nb_modes,
                                                              ((xcFrontDomainLineOpen*)(*it))->getLss3d()),
                                         xApproxFunctionConstant(), key_modifier);
                  break;
               case xcFrontPartBase::LineClose:
                  // space = new xSpaceXFEM(xcSpaceModalFourier(space_name, xSpace::SCALAR, predefined_nb_modes,
                  // ((xcFrontDomainLineClosed*)(*it))->getLss3dCos(), ((xcFrontDomainLineClosed*)(*it))->getLss3dSin()),
                  // xcApproxFunctionLevelSet((*it)->getFr()), key_modifier);
                  space = new xSpaceXFEM(xcSpaceModalFourier(space_name, xSpace::SCALAR, predefined_nb_modes,
                                                             ((xcFrontDomainLineClosed*)(*it))->getLss3dCos(),
                                                             ((xcFrontDomainLineClosed*)(*it))->getLss3dSin()),
                                         xApproxFunctionConstant(), key_modifier);
                  break;
               case xcFrontPartBase::None:
               default:
                  throw;
            }
         }
         break;
         case HERMITE:
         {
            if (minimum_mode_width <= 0.)
            {
               cerr << "Hermite space creation: minimum mode width must be a positive value" << endl;
               throw;
            }
            double front_part_lenght = part->getFrontLenght();
            int nb_modes = (int)(ceil(2 * front_part_lenght / minimum_mode_width));
            if (part->getFrontType() == xcFrontPartBase::LineClose)
            {
               int moduloo = nb_modes % 2;
               if ((moduloo != 0) && (nb_modes > 1)) nb_modes += 1;
            }
            // cout << "xcSpace Hermite nb modes = " << nb_modes << endl;
            (*it)->set_nb_modes(nb_modes);
            switch (part->getFrontType())
            {
               case xcFrontPartBase::Point:
                  // space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR,
                  // (*(part->getFrontRegion().begin(0)))), xcApproxFunctionLevelSet((*it)->getFr()), key_modifier);
                  space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR, (*(part->getFrontRegion().begin(0)))),
                                         xApproxFunctionConstant(), key_modifier);
                  break;
               case xcFrontPartBase::LineOpen:
               {
                  if (nb_modes == 1)
                  {
                     space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR, (*(part->getFrontRegion().begin(0))),
                                                                 &(((xcFrontDomainLineOpen*)(*it))->getLss3d())),
                                            xApproxFunctionConstant(), key_modifier);
                  }
                  else
                  {
                     double mode_width = 4. / (nb_modes - 1);
                     // space = new xSpaceXFEM(xcSpaceModalPolynomeHermiteUniform(space_name, xSpace::SCALAR, nb_modes,
                     // ((xcFrontDomainLineOpen*)(*it))->getLss3d(), mode_width), xcApproxFunctionLevelSet((*it)->getFr()),
                     // key_modifier);
                     space = new xSpaceXFEM(
                         xcSpaceModalPolynomeHermiteUniform(space_name, xSpace::SCALAR, nb_modes,
                                                            ((xcFrontDomainLineOpen*)(*it))->getLss3d(), mode_width),
                         xApproxFunctionConstant(), key_modifier);
                  }
               }
               break;
               case xcFrontPartBase::LineClose:
               {
                  if (nb_modes == 1)
                  {
                     space = new xSpaceXFEM(xcSpaceModalConstant(space_name, xSpace::SCALAR, (*(part->getFrontRegion().begin(0))),
                                                                 &(((xcFrontDomainLineClosed*)(*it))->getLss3dCos())),
                                            xApproxFunctionConstant(), key_modifier);
                  }
                  else
                  {
                     double mode_width = 4. / nb_modes;
                     // space = new xSpaceXFEM(xcSpaceModalPolynomeHermiteUniform(space_name, xSpace::SCALAR, nb_modes,
                     // ((xcFrontDomainLineClosed*)(*it))->getLss3dCos(), ((xcFrontDomainLineClosed*)(*it))->getLss3dSin(),
                     // mode_width), xcApproxFunctionLevelSet((*it)->getFr()), key_modifier);
                     space = new xSpaceXFEM(
                         xcSpaceModalPolynomeHermiteUniform(space_name, xSpace::SCALAR, nb_modes,
                                                            ((xcFrontDomainLineClosed*)(*it))->getLss3dCos(),
                                                            ((xcFrontDomainLineClosed*)(*it))->getLss3dSin(), mode_width),
                         xApproxFunctionConstant(), key_modifier);
                  }
               }
               break;
               case xcFrontPartBase::None:
               default:
                  throw;
            }
            break;
         }
         default:
            cerr << "Not a valid discretisation type" << endl;
            throw;
      }

      spaces.insert(*space);
      delete space;
      ++it;
   }
};

//-------------------------------------------------------------------------------------------------------------------

void xcFrontSpaceManager::exportModes(const std::string& basefilename, int iter)
{
   xSubMesh& globalSubMesh = dmanager->createGlobalSubMesh();
   xRegion region_for_int(&globalSubMesh);
   int dim_region = region_for_int.dim();
   xIntegrationRuleBasic intr_basic;
   xSpace::femFcts fcts;
   xSpace::femKeys keys;
   getSpace().getKeysAndFcts(nullptr, &keys, &fcts);
   xSpace::femFcts::iterator it = fcts.begin();
   xSpace::femFcts::iterator itend = fcts.end();

   xexport::xExportGmshAscii pexport_ascii;
   char filenamemodes[256];
   sprintf(filenamemodes, "Domain_iter_%d", iter);
   xEvalConstant<double> eval_un(1.);
   xexport::Export(eval_un, pexport_ascii, filenamemodes, intr_basic, region_for_int.begin(dim_region),
                   region_for_int.end(dim_region));

   int ifront = 0;
   while (it != itend)
   {
      xEvalApproxFunction<double> eval_modes(*(*it));
      xEvalGradApproxFunction<xtensor::xVector<>> eval_grad_modes(*(*it));

      sprintf(filenamemodes, "%s_part_%d_mod_xx_d_iter_%d", basefilename.c_str(), ifront, iter);
      xexport::Export(eval_modes, pexport_ascii, filenamemodes, intr_basic, region_for_int.begin(dim_region),
                      region_for_int.end(dim_region));
      sprintf(filenamemodes, "%s_grad_part_%d_mod_xx_d_iter_%d", basefilename.c_str(), ifront, iter);
      xexport::Export(eval_grad_modes, pexport_ascii, filenamemodes, intr_basic, region_for_int.begin(dim_region),
                      region_for_int.end(dim_region));
      ifront++;
      it++;
   }
};

//-------------------------------------------------------------------------------------------------------------------

std::vector<std::shared_ptr<xcApproxFunctionLevelSet>> xcFrontSpaceManager::getFrFunctions()
{
   xcFrontDomainManager::iterator it = dmanager->begin_domains_iter();
   xcFrontDomainManager::iterator iten = dmanager->end_domains_iter();
   std::vector<std::shared_ptr<xcApproxFunctionLevelSet>> fr_functions;
   while (it != iten)
   {
      fr_functions.push_back(std::shared_ptr<xcApproxFunctionLevelSet>(new xcApproxFunctionLevelSet((*it)->getFr())));
      it++;
   }
   return fr_functions;
};

//-------------------------------------------------------------------------------------------------------------------

void xcApproxFunctionLevelSet::getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double& res) const
{
   mEntity* e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = ls.getVal(geo_appro->getEntity(), geo_appro->getUVW());
   else
      res = 0.;
}

void xcApproxFunctionLevelSet::getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ,
                                       xtensor::xVector<>& res) const
{
   mEntity* e = geo_appro->getEntity();
   if (ls.getSupport().IsInRegion(e))
      res = ls.getGrad(geo_appro->getEntity());
   else
      res = xtensor::xVector<>(0., 0., 0.);
}
