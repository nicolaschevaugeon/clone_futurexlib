/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// std
#include <cassert>
// xcut
#include "xSupportComponent.h"

using namespace std;

namespace xfem
{
size_t xSupportComponent::getnbComponents() const { return allcomponent.size(); }

const xSupportComponent::component &xSupportComponent::operator()(size_t i) const
{
   assert(i < allcomponent.size());
   return allcomponent[i];
}

xSupportComponent::component &xSupportComponent::operator()(size_t i)
{
   assert(i < allcomponent.size());
   return allcomponent[i];
}

xSupportComponent::component &xSupportComponent::getNewComponent(size_t i)
{
   if (!(i < allcomponent.size()))
   {
      allcomponent.resize(i + 1);
   }
   return allcomponent[i];
}

const xSupportComponent::component &xSupportComponent::getSplitZone() const { return split_zone; }

xSupportComponent::component &xSupportComponent::getSplitZone() { return split_zone; }

void xSupportComponent::addComponent(const component &comp)
{
   allcomponent.push_back(comp);
   return;
}
void xSupportComponent::addComponentByMove(component &comp)
{
   allcomponent.push_back(std::move(comp));
   return;
}

void xSupportComponent::clear()
{
   allcomponent.clear();
   split_zone.clear();
   return;
}

}  // namespace xfem
