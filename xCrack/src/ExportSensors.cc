/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "ExportSensors.h"

#include "xMesh.h"

using AOMD::mEntity;
using xtensor::xPoint;

void xeExportSensors::setMesh(xMesh* mesh)
{
   std::map<int, xfem::xGeomElem*>::iterator itstart = sensor_to_elt.begin();
   //  std::map<int, xfem::xGeomElem*>::iterator itend = sensor_to_elt.end();
   //  sensor_to_elt.erase(itstart,itend);
   for (; itstart != sensor_to_elt.end(); itstart++)
   {
      delete itstart->second;
   }
   sensor_to_elt.clear();

   int i = 0;
   for (const auto& id_pt : sensors)
   {
      xtensor::xPoint p = id_pt.second;
      auto elt_uvw_list = mesh->locateElement(p);
      if (elt_uvw_list.empty())
      {
         cerr << "the sensor located at " << p << " is not in the mesh " << endl;
         abort();
      }
      bool check = false;
      for (auto elt_uvw : elt_uvw_list)
      {
         AOMD::mEntity* pe = const_cast<AOMD::mEntity*>(elt_uvw.first);
         if (test(pe))
         {
            xfem::xGeomElem* geo = new xfem::xGeomElem(pe);
            geo->setUVW(elt_uvw.second);
            if (debug) cout << " in SetMesh, location of sensor is " << geo->getXYZ() << endl;
            sensor_to_elt.insert(make_pair(i, geo));
            check = true;
         }
      }
      if (!check)
      {
         cerr << "the sensor located at " << p << " is is in the mesh but not in a valid element" << endl;
         abort();
      }
      ++i;
   }
}

void xeExportSensors::set(const std::map<int, xtensor::xPoint>& s) { sensors = s; }

void xeExportSensors::read(const xField<>& disp_l, const double& t)
{
   std::map<int, xtensor::xPoint>::const_iterator it = sensors.begin();
   std::map<int, xfem::xGeomElem*>::const_iterator ite = sensor_to_elt.begin();
   for (; it != sensors.end(); ++it, ++ite)
   {
      int is = it->first;
      // xtensor::xPoint p = it->second; -Wunused-but-set-variable
      // xfem::xGeomElem geo(ite->second);
      xfem::xGeomElem* geo = ite->second;
      if (debug) cout << " location of sensor is " << geo->getXYZ() << endl;
      xtensor::xVector<> res;
      disp_l.getVal(geo, geo, res);
      results_disp[is].push_back(make_pair(t, res));
   }
}

/*void xeExportSensors::(const xField<>& disp_l, double &a)
{
  std::map<int, xtensor::xPoint>::const_iterator it = sensors.begin();
  std::map<int, xfem::xGeomElem*>::const_iterator ite = sensor_to_elt.begin();
        std::vector<double> dist;
  for (; it != sensors.end(); ++it, ++ite) {
    xtensor::xPoint p = it->second;
    xfem::xGeomElem* geo = ite->second;
    if (debug) cout << " location of sensor is " << geo->getXYZ() << endl;
    xtensor::xVector<> res;
    disp_l.getVal(geo, geo, res);
                dist.push_back(res[0]);
  }
        a = dist[1]-dist[0];
}*/

void xeExportSensors::measureSensor(const xField<>& field, const int dim, const int ipoint, double& value)
{
   std::map<int, xtensor::xPoint>::const_iterator it = sensors.begin();
   std::map<int, xfem::xGeomElem*>::const_iterator ite = sensor_to_elt.begin();
   int counter = 1;
   bool found = false;

   // the xtensor::xPoint should be given to search in the map...
   for (; it != sensors.end(); ++it, ++ite, counter++)
   {
      if (counter == ipoint)
      {
         found = true;
         break;
      }
   }
   if (found == false) cerr << "ExportSensors: Point " << ipoint << " not found in sensors " << endl;
   // xtensor::xPoint p = it->second; -Wunused-but-set-variable
   xfem::xGeomElem* geo = ite->second;
   if (debug) cout << " location of sensor is " << geo->getXYZ() << endl;
   xtensor::xVector<> res;
   field.getVal(geo, geo, res);
   value = res[dim];
}

void xeExportSensors::print() const
{
   if (debug) cout << " size results of sensors " << results_disp.size() << endl;
   std::map<int, history_t>::const_iterator itd = results_disp.begin();
   std::map<int, xtensor::xPoint>::const_iterator itsensors = sensors.begin();
   for (; itd != results_disp.end(); ++itd, ++itsensors)
   {
      int isensor = itsensors->first;
      std::string filename = "sensors" + int2string(isensor) + ".txt";
      std::ofstream out(filename.c_str());
      const history_t& disp = itd->second;
      history_t::const_iterator itdv = disp.begin();
      for (; itdv != disp.end(); ++itdv)
      {
         double t = itdv->first;
         const xtensor::xVector<>& d = itdv->second;
         out << t << " " << d * t << endl;
      }
      out.close();
   }
}

std::string xeExportSensors::int2string(int i) const
{
   char str[30];
   sprintf(str, "%d", i);
   return std::string(str);
}
