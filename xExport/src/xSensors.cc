/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include <list>
#include <set>

// xfem
#include "xGeomElem.h"
#include "xMesh.h"

// xexport
#include "xSensors.h"

// AOMD
#include "mEntity.h"

// xtool
#include "workInProgress.h"

namespace xexport
{
void xSensors::registerSensor(const int nb_, const xtensor::xPoint pt_, xfem::xMesh* mesh_)
{
   point_container.insert(std::make_pair(nb_, pt_));
   auto list_pe_uvw = mesh_->locateElement(pt_);
   // in // if all proc do have a elts.size()==0 then sensor is out of mesh otherwise it is just in some proc that are not
   // this one. => reduce on max of elts.size to test
   // if (max==0) => message
   // else if (elts.size()==0) => nothing
   // else treatement
   assert(xtool::workInProgress());
   if (list_pe_uvw.empty())
      std::cout << "sensor is out of mesh, it is not set" << std::endl;
   else
   {
      auto pe_uvw = list_pe_uvw.front();
      std::unique_ptr<xfem::xGeomElem> geo(new xfem::xGeomElem(const_cast<AOMD::mEntity*>(pe_uvw.first)));
      geo->setUVW(pe_uvw.second);
      geo->SetIntegrationPointNumberForDegree(0);
      geom_elem_container.insert(std::make_pair(nb_, std::move(geo)));
   }
}

void xSensors::reinit(xfem::xMesh* mesh_)
{
   clean();
   std::vector<int> to_remove;
   for (const auto& id_point : point_container)
   {
      int nb = id_point.first;
      const xtensor::xPoint& p = id_point.second;
      auto list_pe_uvw = mesh_->locateElement(p);
      if (list_pe_uvw.empty())
      {
         std::cout << "sensor " << nb << " is out of mesh, it is removed" << std::endl;
         to_remove.push_back(nb);
      }
      else
      {
         auto pe_uvw = list_pe_uvw.front();
         std::unique_ptr<xfem::xGeomElem> geo(new xfem::xGeomElem(const_cast<AOMD::mEntity*>(pe_uvw.first)));
         geo->setUVW(pe_uvw.second);
         geo->SetIntegrationPointNumberForDegree(0);
         geom_elem_container.insert(std::make_pair(nb, std::move(geo)));
      }
   }
   for (auto id : to_remove) point_container.erase(id);
}

xExportSensors::~xExportSensors()
{
   for (std::map<std::string, std::map<int, std::fstream*>>::iterator it = to_export_container.begin();
        it != to_export_container.end(); ++it)
   {
      std::map<int, std::fstream*>& filestr_map = it->second;
      for (std::map<int, std::fstream*>::iterator itm = filestr_map.begin(); itm != filestr_map.end(); ++itm)
      {
         if (itm->second)
         {
            itm->second->close();
            delete itm->second;
         }
      }
      filestr_map.clear();
   }
   to_export_container.clear();
}

void xExportSensors::init(const xParseData& sensor_parser_, xfem::xMesh* mesh_, const int step_string_length_,
                          std::string key_label_, std::string point_label_)
{
   step_string_length = step_string_length_;
   std::list<std::string> string_list = sensor_parser_.getListString(key_label_);
   for (auto name : string_list)
   {
      to_export_container.insert(std::pair<std::string, std::map<int, std::fstream*>>(name, std::map<int, std::fstream*>()));
   }
   int counter = 1;
   for (auto pt : sensor_parser_.getListVector(point_label_)) sensors.registerSensor(counter++, pt, mesh_);
}

bool xExportSensors::toExport(const std::string& export_name_)
{
   return to_export_container.find(export_name_) != to_export_container.end();
}
}  // namespace xexport
