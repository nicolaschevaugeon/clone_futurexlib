/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _VECTORFIELD_H
#define _VECTORFIELD_H

#include <unordered_map>

#include "mEntity.h"
#include "mVertex.h"
#include "xRegion.h"
#include "xVector.h"

namespace xtensor
{
template <typename T>
class xTensor2;
}

namespace xfem
{
class xVectorField
{
  public:
   //  ... met le drapeau des gradients faux (constructeur par defaut)
   xVectorField() {}
   //  ... met la levelset a zero et le drapeau des gradients faux
   xVectorField(const xRegion& support, const xtensor::xVector<>& val = xtensor::xVector<>(0., 0., 0.));

   template <class T>
   void load(const T& func)  //; //the T class must define operator ()(const xtensor::xPoint&)
   {
      for (xIter it = support.begin(0); it != support.end(0); ++it)
      {
         AOMD::mVertex* v = (AOMD::mVertex*)*it;
         f[v] = func(v->point());
      }
   }

   void setSupport(const xRegion& m, const xtensor::xVector<>& val = xtensor::xVector<>(0., 0., 0.));

   // set Function
   xtensor::xVector<>& operator()(AOMD::mEntity* e);

   // get Functions
   const xtensor::xVector<>& operator()(AOMD::mEntity* e) const;

   // get vals aux noeuds
   std::vector<xtensor::xVector<>> getVals(AOMD::mEntity* e) const;

   xtensor::xVector<> getVal(AOMD::mEntity* e, const xtensor::xPoint& uvw) const;

   xtensor::xTensor2<double> getGrad(AOMD::mEntity* e) const;

   // get mesh
   xRegion getSupport() const;
   void clear();

   void exportGmsh(const std::string& name);
   void exportGmsh(const std::string& name, xRegion sub);

   void exportMatlab(std::ostream& fout, const std::string& field_name, int level = 0);

   void printDebug() const;

   void importGmsh(const std::string& filename);

  private:
   mutable xRegion support;
   std::unordered_map<AOMD::mEntity*, xtensor::xVector<>, AOMD::EntityHashKey, AOMD::EntityEqualKey> f;
};

}  // namespace xfem

#endif
