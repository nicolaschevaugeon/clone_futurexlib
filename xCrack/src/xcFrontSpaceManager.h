/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xcFrontSpaceManager_
#define _xcFrontSpaceManager_

#include <string>
#include <vector>
// xfem
#include "xApproxFunction.h"
#include "xLevelSet.h"
// xcrack
#include "xcFrontDomainManager.h"
#include "xcFrontSpaces.h"

enum discreteType
{
   CONSTANT,
   MODAL,
   HERMITE
};
// "constant" is the constant mode,
// "modal" is a modal representation with a non zero value of the modes on all ls1d,
// (legendre modes are used for open lines, fourier modes for closed lines)
// "hermite" is a partition of unity representation using hermite polynomials

class xcApproxFunctionLevelSet : public xfem::xApproxFunction
{
  public:
   xcApproxFunctionLevelSet(const xfem::xLevelSet& l) : ls(l) {}
   void getVal(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, double&) const override;
   void getGrad(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<>&) const override;
   std::string name() override { return "xcApproxFunctionLevelSet"; }

  private:
   const xfem::xLevelSet& ls;
};

class xcFrontSpaceManager
{
  public:
   xcFrontSpaceManager(xcFrontDomainManager* _dmanager, int discretisationType, const std::string& space_name,
                       double minimum_mode_width = 0.);
   ~xcFrontSpaceManager();
   xfem::xSpaceComposite& getSpace() { return spaces; }
   std::vector<std::shared_ptr<xcApproxFunctionLevelSet>> getFrFunctions();
   xcFrontDomainManager* dmanager;
   void exportModes(const std::string& basefilename, int iter);

  private:
   xfem::xSpaceComposite spaces;
   //  std::vector<xApproxFunctionLevelSet > fr_functions;
};

#endif
