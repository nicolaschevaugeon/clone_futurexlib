/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#ifndef XFEM_XSUPPORTCOMPONENT_H
#define XFEM_XSUPPORTCOMPONENT_H

// std
#include <functional>
#include <vector>

// AOMD
#include "mAttachableDataContainer.h"
namespace AOMD
{
class mEntity;
}

namespace xcut
{
class xPhysDomain;
}
// xfem

namespace xfem
{
class xSupportComponent
{
  public:
   typedef std::vector<const AOMD::mEntity *> component;

   size_t getnbComponents() const;

   const component &operator()(size_t i) const;
   component &operator()(size_t i);

   const component &getSplitZone() const;
   component &getSplitZone();

   void addComponent(const component &comp);
   void addComponentByMove(component &comp);

   void clear();
   component &getNewComponent(size_t i);

  private:
   std::vector<component> allcomponent;
   component split_zone;
};

// a typedef to discribe functor that return a xSupportComponent from a entity
typedef std::function<const xfem::xSupportComponent &(AOMD::mEntity *)> getSupportComponent_functor_t;

}  // end namespace xfem
#endif
