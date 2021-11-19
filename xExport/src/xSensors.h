/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef _SENSOR_H
#define _SENSOR_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

// xtensor
#include "xPoint.h"
#include "xTensor2.h"

// xfem
#include "xCommandOnGeomElem.h"
#include "xEval.h"
#include "xIntegrationRule.h"
#include "xParseData.h"
#include "xRegion.h"

// A xSensors is a collection of xPoint numbered (map) where measure may
// be needed. Measure is done through the use of xEval<T>.
// REMARKS:
// (1) nothing is done in the constructor, the "set" function add new xPoint
// (2) xPoint is checked whether or not it is on the mesh,
// if not it is only removed (no assert())
// (3) "measure" measures only on xPoint number nb_ (and check if it exists, if not return false)
// (4) "measureAll" measures all existing xPoint

// A xExportSensors has a member xSensors and works with the same behavior.
// REMARKS:
// (1) it writes values into files "on the fly".
// (2) files names must be defined in a .dat file read from xParseData EXPORT_SENSOR_LABELS, EXPORT_SENSOR_POINTS.
// (3) if a label in "measure" or "measureAll" doesn't match an EXPORT_SENSOR_LABELS, then no measure is done.

namespace xexport
{
class xSensors
{
  public:
   xSensors() {}
   ~xSensors() { clean(); }

   // set new mPoint to the point_container, and xfem::xGeomElem* to the geom_elem_container
   void registerSensor(const int nb_, const xtensor::xPoint pt_, xfem::xMesh* mesh_);

   void registerBox(const int nb_, const xfem::xEntityFilter box_) { box_filter_container.insert(std::make_pair(nb_, box_)); }

   // remove mPoint from the point_container and xfem::xGeomElem* from the geom_elem_container
   inline void removeSensor(const int nb_);

   // remove all point from the point_container and all xfem::xGeomElem* from the geom_elem_container
   inline void removeAll();

   // reinit all xfem::xGeomElem* of the geom_elem_container from point_container
   void reinit(xfem::xMesh* mesh_);

   // measure value from an evaluator on mPoint number nb_
   template <typename T>
   bool measure(const xfem::xEval<T>& eval_, const int nb_, T& val_);

   // measure value from an evaluator on all mPoint
   template <typename T>
   void measureAll(const xfem::xEval<T>& eval_, std::map<int, T>& val_map_);

   // measure maximum value on a predefined box
   template <typename ITER>
   bool measureMaxOnBox(const xfem::xEval<double>& eval_, const int nb_, const xfem::xIntegrationRule& integration_rule_,
                        const ITER& begin_, const ITER& end_, double& val_)
   {
      std::map<int, xfem::xEntityFilter>::const_iterator found = box_filter_container.find(nb_);
      if (found != box_filter_container.end())
      {
         val_ = -std::numeric_limits<double>::max();
         xfem::xEvalMaxCommand<const xfem::xEval<double>> command(eval_, val_);
         xfem::xFilteredRegion<ITER, xfem::xEntityFilter> fr(begin_, end_, found->second);
         ApplyCommandOnIntegrationRule(command, integration_rule_, fr.begin(), fr.end());
         return true;
      }
      return false;
   }

   template <typename ITER>
   void measureMaxOnAllBoxes(const xfem::xEval<double>& eval_, const xfem::xIntegrationRule& integration_rule_,
                             const ITER& begin_, const ITER& end_, std::map<int, double>& val_map_)
   {
      double val = 0.;
      for (const auto& id_filter : box_filter_container)
      {
         int id = id_filter.first;
         measureMaxOnBox(eval_, id, integration_rule_, begin_, end_, val);
         if (val > -std::numeric_limits<double>::max())
         {
            val_map_.insert(std::pair<int, double>(id, val));
         }
      }
   }

  private:
   // clean geom_elem_container
   inline void clean();

   std::map<int, xtensor::xPoint> point_container;
   std::map<int, std::unique_ptr<xfem::xGeomElem>> geom_elem_container;
   std::map<int, xfem::xEntityFilter> box_filter_container;
};

class xExportSensors
{
  public:
   xExportSensors() { xtensor::xTensor2<>::setExportFormat(0); }
   ~xExportSensors();

   // read xParseData EXPORT_SENSOR_LABELS to get all file names and EXPORT_SENSOR_POINTS to get all xPoint
   void init(const xParseData& sensor_parser_, xfem::xMesh* mesh_, const int step_string_length_ = 2,
             std::string key_label_ = "EXPORT_SENSOR_LABELS", std::string point_label_ = "EXPORT_SENSOR_POINTS");

   bool toExport(const std::string& export_name_);

   // reinits xSensor geom_elem_container member
   inline void reinit(xfem::xMesh* mesh_);

   // measures and writes into a file named
   // "sensor_"+nb_+export_name_+".txt" if sensor nb_ exists and
   // if export_name_ matches a string in to_export_container
   template <typename T>
   void measure(const xfem::xEval<T>& eval_, const double time_, const int nb_, const std::string export_name_);

   // measures and writes into files named
   // "sensor_"+nb+export_name+".txt" for all existing sensor nb
   // if export_name_ matches a string in to_export_container
   template <typename T>
   void measureAll(const xfem::xEval<T>& eval_, const double time_, const std::string export_name_);

   template <typename T>
   void measure(const T& val_, const double time_, const std::string export_name_);

  private:
   inline std::string outputId(const int time_id);

   xSensors sensors;
   std::map<std::string, std::map<int, std::fstream*>> to_export_container;
   int step_string_length;
};

void xSensors::removeSensor(const int nb_)
{
   geom_elem_container.erase(nb_);
   point_container.erase(nb_);
}

void xSensors::removeAll()
{
   clean();
   point_container.clear();
}

void xSensors::clean() { geom_elem_container.clear(); }

template <typename T>
bool xSensors::measure(const xfem::xEval<T>& eval_, const int nb_, T& val_)
{
   auto found = geom_elem_container.find(nb_);
   if (found != geom_elem_container.end())
   {
      const xfem::xGeomElem* geo = found->second.get();
      eval_(geo, geo, val_);
      return true;
   }
   else
   {
      std::cout << "sensor " << nb_ << " not defined or not in the mesh" << std::endl;
      return false;
   }
}

template <typename T>
void xSensors::measureAll(const xfem::xEval<T>& eval_, std::map<int, T>& val_map_)
{
   for (const auto& id_geom_elem : geom_elem_container)
   {
      int nb = id_geom_elem.first;
      xfem::xGeomElem* geo = id_geom_elem.second.get();
      T val;
      eval_(geo, geo, val);
      val_map_.insert(std::pair<int, T>(nb, val));
   }
}

std::string xExportSensors::outputId(const int time_id)
{
   std::ostringstream out;
   out.width(step_string_length);
   out.fill('0');
   out << time_id;
   return "_" + out.str();
}

void xExportSensors::reinit(xfem::xMesh* mesh_) { sensors.reinit(mesh_); }

template <typename T>
void xExportSensors::measure(const xfem::xEval<T>& eval_, const double time_, const int nb_, const std::string export_name_)
{
   std::map<std::string, std::map<int, std::fstream*>>::iterator found = to_export_container.find(export_name_);
   if (found != to_export_container.end())
   {
      T val;
      if (sensors.measure(eval_, nb_, val))
      {
         std::map<int, std::fstream*>& filestr_map = found->second;
         std::map<int, std::fstream*>::iterator found_filestr = filestr_map.find(nb_);
         if (found_filestr != filestr_map.end())
         {
            *(found_filestr->second) << time_ << " " << val << std::endl;
         }
         else
         {
            std::string file = "sensor" + outputId(nb_) + "_" + export_name_ + ".txt";
            std::fstream* filestr = new std::fstream(file.c_str(), std::fstream::out);
            filestr_map.insert(std::pair<int, std::fstream*>(nb_, filestr));
            *filestr << time_ << " " << val << std::endl;
         }
      }
   }
}

template <typename T>
void xExportSensors::measureAll(const xfem::xEval<T>& eval_, const double time_, const std::string export_name_)
{
   std::map<std::string, std::map<int, std::fstream*>>::iterator found = to_export_container.find(export_name_);
   if (found != to_export_container.end())
   {
      std::map<int, T> val_map;
      sensors.measureAll(eval_, val_map);

      std::map<int, std::fstream*>& filestr_map = found->second;

      for (const auto& id_val : val_map)
      {
         int nb = id_val.first;
         T val = id_val.second;
         std::map<int, std::fstream*>::iterator found_filestr = filestr_map.find(nb);
         if (found_filestr != filestr_map.end())
         {
            *(found_filestr->second) << time_ << " " << val << std::endl;
         }
         else
         {
            std::string file = "sensor" + outputId(nb) + "_" + export_name_ + ".txt";
            std::fstream* filestr = new std::fstream(file.c_str(), std::fstream::out);
            filestr_map.insert(std::pair<int, std::fstream*>(nb, filestr));
            *filestr << time_ << " " << val << std::endl;
         }
      }
   }
}

template <typename T>
void xExportSensors::measure(const T& val_, const double time_, const std::string export_name_)
{
   std::map<std::string, std::map<int, std::fstream*>>::iterator found = to_export_container.find(export_name_);
   if (found != to_export_container.end())
   {
      std::map<int, std::fstream*>& filestr_map = found->second;
      std::map<int, std::fstream*>::iterator found_filestr = filestr_map.find(0);
      if (found_filestr != filestr_map.end())
      {
         *(found_filestr->second) << time_ << " " << val_ << std::endl;
      }
      else
      {
         std::string file = "sensor_val_" + export_name_ + ".txt";
         std::fstream* filestr = new std::fstream(file.c_str(), std::fstream::out);
         filestr_map.insert(std::pair<int, std::fstream*>(0, filestr));
         *filestr << time_ << " " << val_ << std::endl;
      }
   }
}
}  // namespace xexport

#endif
