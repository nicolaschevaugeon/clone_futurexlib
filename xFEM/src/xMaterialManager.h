/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __material__manager_h
#define __material__manager_h

#include <string>
#include "xTensorsPtr.h"
#include "xSingleton.h"
#include "xFactory.h"
#include "xMaterial.h"
#include "xZone.h"
#include "xGeomElem.h"
#include "xTensors.h"
#include "xVariabManager.h"	


namespace xfem
{

class xMaterialManager;
typedef xtool::xSingleton<xMaterialManager> xMaterialManagerSingleton;


template <class T> 
class xMaterialCreator
{ 
public: 
  xMaterial* operator()() const { return new T; } //VALGRINDLEAK 
};


class  xTensorsValueCreator;
 
class xMaterialManager 
{

public:

xMaterialManager();
~xMaterialManager();
zonePtr createZone(const std::string& zone_name, const std::string& mat_class, const std::string& mat_param);

zonePtr createZone(int entity, const std::string& mat_class, const std::string& mat_param);

xMaterial* createMaterial(const std::string& mat_class, const std::string& mat_param);

xMaterial* getMaterial(const xGeomElem* geo_integ);

xMaterial* getMaterialForZone(int zone_id);

xMaterial* getMaterialForZone(const std::string& zone_name);

xValue< tensorsPtr_t >* getMaterialVariables(const xGeomElem* geo_integ, xVariabManager& v);

xValue< tensorsPtr_t >* createMaterialVariables(const xGeomElem* geo_integ, xVariabManager& v,
					       xTensorsValueCreator&  vc);

void deleteMaterialVariables(const xGeomElem* geo_integ, xVariabManager& v);

bool registerMaterial(const std::string& id, std::function<xMaterial* (void)> creator)
{
  return factory.registerCreator(id, creator);
}

private:
 xMaterial*       curr_mat;
 int              curr_zone;
 xZoneContainer*  zones;
 std::vector<xValKey::ids_size_t> integration_point_id_map;
 xValKey::ids_size_t get_integration_point_id_map( size_t idgauss);
 xFactory<xMaterial, std::string, std::function<xMaterial* (void)> > factory;
 std::list<xMaterial *> listcreatedmat;
};

} // end of namespace

#endif
