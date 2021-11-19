/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xZone.h"

#include <iostream>
#include <string>

#include "xMaterial.h"
//

namespace xfem
{
using std::cerr;
using std::endl;
using std::string;

xtool::xStringManager<int> xZone::names;

xZone::xZone() : id(-1), zone_name(""), material(nullptr) {}
xZone::~xZone() = default;
xZone::xZone(int i, xMaterial* mat) : id(i), zone_name(names.getName(i)), material(mat) {}
xZone::xZone(const string& z, xMaterial* mat) : id(names.getId(z)), zone_name(z), material(mat) {}
zonePtr xZoneContainer::getZoneWithId(int id)
{
   iterator found = zones.find(id);
   if (found == zones.end())
   {
      cerr << "no zone with id " << id << endl;
      assert(0);
   }
   return found->second;
}
void xZoneContainer::addZone(zonePtr zone) { zones.insert(std::make_pair(zone->GetId(), zone)); }

}  // namespace xfem
