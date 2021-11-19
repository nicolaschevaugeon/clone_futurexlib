/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _ExportSensors_hhh_
#define _ExportSensors_hhh_

#include <complex>

#include "xData.h"
#include "xField.h"
#include "xPointToDouble.h"

using namespace xfem;

class xeTrue
{
  public:
   bool operator()(AOMD::mEntity* e) const { return true; }
};

// PEB: needs to be more general and probably moved to Xfem...!
class xeExportSensors
{
  public:
   xeExportSensors(int p_, const std::map<int, xtensor::xPoint>& points,
                   const std::function<bool(AOMD::mEntity* e)>& t_ = xeTrue(), bool _deb = false)
       : period(p_), test(t_), debug(_deb)
   {
      set(points);
   }
   bool exportRequired(int id) const
   {
      if (debug) std::cout << " the test for export sensors required is " << (id % period == 0) << std::endl;
      return (id % period == 0);
   }
   ~xeExportSensors()
   {
      for (auto& id_geom : sensor_to_elt) delete id_geom.second;
   }

   void setMesh(xMesh* m);
   void read(const xField<>&, const double& t);
   void print() const;
   inline int size() const { return sensors.size(); }

   // returns the value of the dim_th component of the field at ipoint_th point
   void measureSensor(const xField<>& field, const int dim, const int ipoint, double& value);

   typedef std::vector<std::pair<double, xtensor::xVector<>>> history_t;

   template <class TYPE>
   void measureSensor(const xEval<TYPE>& eval, const int dim, const int ipoint, TYPE& value)
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
      if (found == false) std::cerr << "ExportSensors: Point " << ipoint << " not found in sensors " << std::endl;
      xtensor::xPoint p = it->second;
      xfem::xGeomElem* geo = ite->second;
      if (debug) std::cout << " location of sensor is " << geo->getXYZ() << std::endl;
      eval(geo, geo, value);
   }

  private:
   void set(const std::map<int, xtensor::xPoint>& s);
   int period;
   std::map<int, xtensor::xPoint> sensors;
   std::map<int, xfem::xGeomElem*> sensor_to_elt;
   // first int is for sensor id, second for time
   std::map<int, history_t> results_disp;
   std::string int2string(int i) const;
   std::function<bool(AOMD::mEntity* e)> test;
   bool debug;
};

#endif
