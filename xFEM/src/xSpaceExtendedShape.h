/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef XSPACE_EXTEND_SHAPE_H
#define XSPACE_EXTEND_SHAPE_H
#include <map>
#include <vector>

#include "xSpace.h"

namespace xfem
{
class xExtendShapeGeneratorBase;
class xExtendShapeFcts;

class xSpaceExtendedShape : public xSpaceRegular
{
  public:
   ~xSpaceExtendedShape() override;
   xSpaceExtendedShape(const std::string& space_name, const std::string& physical_name,
                       xfem::xExtendShapeGeneratorBase* generator_);
   void getKeys(AOMD::mEntity* e, femKeys* keys) override;
   void getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* appro) override;

  private:
   xfem::xExtendShapeGeneratorBase* generator;
   const xfem::xExtendShapeFcts* curent_extend;
   void setCurentExtend(AOMD::mEntity* e);
};

}  // namespace xfem

#endif
