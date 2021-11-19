/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __SPACEFACTORY_H
#define __SPACEFACTORY_H

#include <string>

#include "xSpace.h"
#include "xSpacePolynomial.h"
#include "xSpacePolynomialQH.h"

namespace xfem
{
class xNonLocalInfoGeneratorForKeysAndFcts;

class xSpaceFactoryBase
{
  public:
   enum space_product_tensorial_type
   {
      SCALAR,
      V1Dx,
      V1Dy,
      V1Dz,
      V2Dxy,
      V2Dxz,
      V2Dyz,
      V3D
   };

   void setSpaceProductOrder(const int order_);
   void setSpaceProductPhysStrings(std::vector<std::string> &phys_strings_);
   void setSpaceProductGenerator(xNonLocalInfoGeneratorForKeysAndFcts *generator);
   const std::string &getPhysString(unsigned int id) const;
   virtual ~xSpaceFactoryBase() = default;
   virtual xSpace::spacePtr getSpace() = 0;
   virtual xSpace::spacePtr getSpace(std::string s) = 0;

  protected:
   template <typename SPACE>
   SPACE createSpaceProduct(std::string &phys, xSpace::TensorialType_t space_tensorial);

   xSpaceLagrange::lag_degree_t getHierachicalOrder();

   void setPhysString(std::string s1);
   void forceSetPhysString(std::string &s1);
   void setPhysString(std::string s1, std::string s2);
   void setPhysString(std::string s1, std::string s2, std::string s3);

   void error1();

   xNonLocalInfoGeneratorForKeysAndFcts *generator;
   std::vector<std::string> phys_strings;
   int order;
};
template <typename SPACE, xSpaceFactoryBase::space_product_tensorial_type TT>
class xSpaceFactory : public xSpaceFactoryBase
{
  public:
   xSpaceFactory();
   xSpace::spacePtr getSpace() override;
   xSpace::spacePtr getSpace(std::string s) override;

  private:
};

}  // namespace xfem

#include "xSpaceFactory_imp.h"

#endif
