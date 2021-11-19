/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <string>
#include <cstdio>
#include <iostream>
#include "xZone.h"
#include "xValKey.h"
#include "xGeomElem.h"
#include "xMaterialManager.h"
#include "xValue.h"
#include "xValueCreators.h"
using namespace std;

namespace xfem
{

xMaterialManager::xMaterialManager() 
: curr_mat(nullptr), curr_zone(-1), zones(new xZoneContainer) 
{ 
  /*
  for (int i = 0; i < 2000; ++i) {
   char str[10];
   std::sprintf(str, "%d", i+1);
   std::string id(str);
   std::string name("INTEG_POINT_" + id);
   integration_point_id_map.push_back(xKeyInfo::getGeomId(name)); 
 }
  */
}

  xMaterialManager::~xMaterialManager(){
    // std::cout << listcreatedmat.size() << std::endl;
    for ( std::list<xMaterial *>::iterator it = listcreatedmat.begin(); it != listcreatedmat.end(); ++it){
      //    std::cout << *it << std::endl;
       delete *it;
    }
    listcreatedmat.clear();
    delete zones;
  }

xValKey::ids_size_t xMaterialManager::get_integration_point_id_map( size_t idgauss){
  if  (integration_point_id_map.size()== 0){
    integration_point_id_map.reserve(2000);
    for (int i = 0; i < 2000; ++i) {
      char str[10];
      std::sprintf(str, "%d", i+1);
      std::string id(str);
      std::string name("INTEG_POINT_" + id);
      integration_point_id_map.push_back(xKeyInfo::getGeomId(name)); 
    }
  }
  if (idgauss < 2000) return integration_point_id_map[idgauss];
  else throw;
}

zonePtr xMaterialManager::createZone(int zone_entity, 
				      const string& mat_class, 
				      const string& filename) 
{
  const bool debug = xdebug_flag;
  if (debug) cout << "mat_class in create_zone " << mat_class 
                  << " filename in create_zone " << filename << endl;
  zonePtr  zone_ptr = zonePtr(new xZone(zone_entity, createMaterial(mat_class, filename)));
  zones->addZone(zone_ptr);
  return zone_ptr;
}

zonePtr xMaterialManager::createZone(const string& zone_name, 
					const string& mat_class, 
					const string& filename) 
{
  zonePtr  zone_ptr = zonePtr(new xZone(zone_name, createMaterial(mat_class, filename)));
  zones->addZone(zone_ptr);
  return zone_ptr;
}


xMaterial* xMaterialManager::createMaterial(const string& mat_class, const string& filename) 
{  
  const bool debug =  xdebug_flag;
  if (debug) cout << " creating material " << mat_class << " from file " <<  filename << endl;
  xMaterial* mat =  factory.createObject(mat_class);
  if (debug) cout << " material created " << mat << std::endl;
  mat->readProperties(filename);
  mat->checkProperties();
  listcreatedmat.push_back(mat);
  return mat;
}


xMaterial* xMaterialManager::getMaterial(const xGeomElem* geo_integ)
{ 
  const bool debug = xdebug_flag;
  if (debug) cout << " in getMaterial(GeomElem_c* geo_integ) " << endl;
  if (debug) cout << curr_zone << " curr_zone " << endl;
  if (curr_zone != geo_integ->getZone()){
    if (debug) cout << "before getZone " << curr_zone;
    curr_zone = geo_integ->getZone();
    if (debug) cout << "zone_id is " << curr_zone;
    zonePtr zone = zones->getZoneWithId(curr_zone);
    if (debug) assert(zone);
    curr_mat = zone->GetMaterial();
    if (debug) assert(curr_mat);
  }
  return curr_mat;
}

xMaterial* xMaterialManager::getMaterialForZone(int zone_id)
{ 
  const bool debug = xdebug_flag;
  if (curr_zone != zone_id){
    curr_zone = zone_id;
    if (debug) cout << "zone_id is " << curr_zone;
    zonePtr zone = zones->getZoneWithId(curr_zone);
    if (debug) assert(zone);
    curr_mat = zone->GetMaterial();
    if (debug) assert(curr_mat);
  }
  return curr_mat;
}


xMaterial* xMaterialManager::getMaterialForZone(const string& zone_name)
{ 
  //const bool debug = xdebug_flag;
  int zone_id = xZone::names.getId(zone_name);
  return getMaterialForZone(zone_id);
}

xValue< tensorsPtr_t >* xMaterialManager::getMaterialVariables(const xGeomElem* geo_integ, 
									      xVariabManager& var_manager)
{
  xMaterial* mat = getMaterial(geo_integ);
  const xTensors* properties = mat->getProperties();
  // xValKey key(xKeyInfo::getPhysId(properties->astring("NAME")),
  //	       integration_point_id_map[geo_integ->getCurrentIntegrationPointId()],
  //	       geo_integ->getEntity());     

  xValKey key(xKeyInfo::getPhysId(properties->astring("NAME")),
	      get_integration_point_id_map(geo_integ->getCurrentIntegrationPointId()),
	      geo_integ->getEntity()); 
  
  return var_manager.find(key);
}

xValue< tensorsPtr_t >* xMaterialManager::createMaterialVariables(const xGeomElem* geo_integ, 
										xVariabManager& var_manager,
										xTensorsValueCreator& vc)
{
  
  const bool debug = xdebug_flag;
  xMaterial* mat = getMaterial(geo_integ);
  const xTensors* properties = mat->getProperties();
  vc.setSignature(mat->getVariablesSignature());
  if (debug) cout << "point _id " <<  geo_integ->getCurrentIntegrationPointId() << endl;
  if (debug) cout << " corresponding integer " <<  
	       //   integration_point_id_map[geo_integ->getCurrentIntegrationPointId()] << endl;
	       get_integration_point_id_map(geo_integ->getCurrentIntegrationPointId()) << endl;
  xValKey key(xKeyInfo::getPhysId(properties->astring("NAME")),
	      //integration_point_id_map[geo_integ->getCurrentIntegrationPointId()],
	      get_integration_point_id_map(geo_integ->getCurrentIntegrationPointId()),
	      geo_integ->getEntity());     
  return var_manager.insert(key, vc);
}

void xMaterialManager::deleteMaterialVariables(const xGeomElem* geo_integ, 
                                               xVariabManager& var_manager)
{
  xMaterial* mat = getMaterial(geo_integ);
  const xTensors* properties = mat->getProperties();
  xValKey key(xKeyInfo::getPhysId(properties->astring("NAME")),
	      get_integration_point_id_map(geo_integ->getCurrentIntegrationPointId()),
	      geo_integ->getEntity());
  var_manager.erase(key);
}


} // end of namespace
