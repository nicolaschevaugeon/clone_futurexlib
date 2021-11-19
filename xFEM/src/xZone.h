/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __ZONE_H
#define __ZONE_H

#include <map>
#include <string>
#include <memory>
#include "xStringManager.h"

namespace xfem
{

class xZone;
class xMaterial;

using  zonePtr = std::shared_ptr<xZone>;

class xZone {
private:
  int id;
  std::string zone_name;
  xMaterial* material;
public:
  xZone();
  xZone(int i, xMaterial* mat);
  xZone(const std::string& name, xMaterial* mat);
  ~xZone();  
  int GetId()                     const {return id; }
  xMaterial *GetMaterial()       const {return material;}
  static xtool::xStringManager<int> names;  
};

class xZoneContainer {
private :
   typedef std::map<int, zonePtr> rep_type;
   rep_type zones;
public :   
  typedef rep_type::const_iterator const_iterator;
  typedef rep_type::const_iterator       iterator;
  typedef rep_type::value_type         value_type; 
  zonePtr getZoneWithId(int id);
  void addZone(zonePtr zone);
};


} // end of namespace













#endif






