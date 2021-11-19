/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef _GetJint_H
#define _GetJint_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

// xtensor
#include "xTensor2.h"
#include "xTensor4.h"
#include "xVector.h"

// xfem
#include "xField.h"
#include "xGeomElem.h"
#include "xRegion.h"
#include "xValueManager.h"
#include "xVariabManager.h"

// xcrack
#include "CrackPostpro.h"

class lCrack;

using namespace AOMD;
using namespace xfem;

namespace xcrack
{
class J3DCommand_c : public xCommandOnGeomElem
{
  public:
   J3DCommand_c(const xField<>& q, const xField<>& disp, const lCrack& c, xField<>& j);
   void openApproxElem(xfem::xGeomElem* g_appro) override;
   void closeApproxElem(xfem::xGeomElem* g_appro) override;
   void execute(xfem::xGeomElem* geo_integ) override;
   xTensors getTensors();

  protected:
   xTensorsSignature signature;
   xTensors tensors;
   double alpha;
   xtensor::xVector<> dalpha;
   xtensor::xTensor4<> DC;
   const lCrack& crack;

   std::vector<double> MatInfo;
   xEvalField<xtool::xIdentity<double>> eval_alpha;
   xEvalGradField<xtool::xIdentity<xtensor::xVector<>>> grad_alpha;
   xEvalGradField<xtool::xIdentity<xtensor::xTensor2<>>> grad_disp;
   xField<>& j_info;
};

template <template <class> class DATAMANAGER = xMesh::datamanager_t>
class xFind3DAroundFront
{
  public:
   xFind3DAroundFront() : was_created_by(xMesh::get_const_was_created_by()) {}

   AOMD::mEntity* operator()(AOMD::mEntity* e)
   {
      mEntity* eup = was_created_by.at(*e);
      mEntity* e3d = was_created_by.at(*eup);
      if (e3d->getLevel() != 3) return e3d->get(3, 0);
      return e3d;
   }

  private:
   const DATAMANAGER<AOMD::mEntity*>& was_created_by;
};

class xAcceptFrontInBox
{
  public:
   xAcceptFrontInBox(xRegion& b) : box(b) {}
   bool operator()(AOMD::mEntity* e)
   {
      if (box.find(to3d(e))) return true;
      return false;
   }

  private:
   xFind3DAroundFront<> to3d;
   xRegion& box;
};

}  // namespace xcrack
#endif
