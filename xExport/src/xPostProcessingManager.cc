/*
    This file is a part of eXlibris C++ Library
    under the GNU General Public License:
    See the LICENSE.md files for terms and
    conditions.
*/

#include "xPostProcessingManager.h"

#include "mEntity.h"
#include "xFiniteElement.h"

using namespace std;
using namespace xfem;
using namespace xexport;

int GetNbDigits(int nb)
{
   if (nb < 0)
   {
      return 1;
   }
   int digits = 0;
   while (nb)
   {
      nb /= 10;
      ++digits;
   }
   return digits;
}

xPostProcessingManager::xPostProcessingManager(xexport::xExport& pexport_, xEntityFilter filter)
    : filter(filter), pexport(pexport_)
{
}

void xPostProcessingManager::initExportManager(const xParseData& parse_data_, std::string export_manager_label,
                                               const int nb_step_max)
{
   export_manager.init(parse_data_, export_manager_label);
   export_manager.setNumerationSize(GetNbDigits(nb_step_max));
}

void xPostProcessingManager::initSensors(xfem::xMesh& mesh, const xParseData& parse_data_, std::string export_sensors_label,
                                         std::string export_sensors_point)
{
   sensor_manager.init(parse_data_, &mesh, 2, export_sensors_label, export_sensors_point);
}
void xPostProcessingManager::reinitMesh(xMesh& mesh) { sensor_manager.reinit(&mesh); }

void xPostProcessingManager::setDebugFlag(bool debug_) { debug = debug_; }

std::string xPostProcessingManager::getFileName() const { return export_manager.getFileName(); }
bool xPostProcessingManager::toExport(const std::string export_name) { return sensor_manager.toExport(export_name); }

bool xPostProcessingManager::toExport(const std::string export_name, const int step, const std::string& extension_name)
{
   return export_manager.toExport(export_name, step, extension_name);
}

void xPostProcessingManager::exportOnSpace(const std::string export_name, const int step, const xLevelSet& ls)
{
   if (export_manager.toExport(export_name, step, ""))
   {
      if (debug) cout << "Exporting " << export_manager.getFileName() << " with nb split = " << pexport.getNbSplit() << endl;
      pexport.openFile(export_manager.getFileName());
      Export(ls, pexport, export_name);
      pexport.closeFile();
      pexport.setNbSplit(1);
   }
}

void xPostProcessingManager::exportOnSpace(const std::string export_name, const int step, const xMesh& mesh)
{
   if (export_manager.toExport(export_name, step, ""))
   {
      if (debug) cout << "Exporting " << export_manager.getFileName() << " with nb split = " << pexport.getNbSplit() << endl;
      Export(mesh, pexport, export_manager.getFileName());
      pexport.setNbSplit(1);
   }
}

void xPostProcessingManager::exportOnTime(const std::string name, const double time, const double val)
{
   sensor_manager.measure(val, time, name);
}

void xPostProcessingManager::exportOnTime(const std::string name, const double time, const xEval<double>& eval)
{
   sensor_manager.measureAll(eval, time, name);
}

void xPostProcessingManager::exportOnTime(const std::string name, const double time, const xEval<xtensor::xVector<>>& eval)
{
   sensor_manager.measureAll(eval, time, name);
}

void xPostProcessingManager::exportOnTime(const std::string name, const double time, const xEval<xtensor::xTensor2<>>& eval)
{
   sensor_manager.measureAll(eval, time, name);
}

void xPostProcessingManager::exportDebug(const std::string export_name, const int step, xfem::xValueManagerDist<double>& dm)
{
   if (export_manager.toExport(export_name, step, ".dbg"))
   {
      if (debug) cout << "Exporting " << export_manager.getFileName() << endl;
      dm.PrintForDebug(export_manager.getFileName());
   }
}

void xPostProcessingManager::saveMesh(std::string export_name, int step, xMesh& mesh)
{
   if (step > 0 && export_manager.toExport("save", step, ""))
   {
      Export(mesh, pexport, export_name + "_save");
   }
}

void xPostProcessingManager::saveField(std::string export_name, int step, const xField<>& field, xRegion reg)
{
   if (step > 0 && export_manager.toExport("save", step, ""))
   {
      xValueManagerDist<double>* double_manager = field.getValueManager();
      std::ofstream oss(export_name + "_save.bin", ios::binary);
      if (oss.is_open())
      {
         const int info_size = 2 * reg.dim() * sizeof(int);
         int nb = reg.size(reg.dim());
         oss.write((char*)&nb, sizeof(int));
         for (auto e : reg)
         {
            xFiniteElement FEM;
            FEM.setKeys(e, field.begin(), field.end());
            const int size = FEM.sizeKey();
            std::vector<xValue<double>*> vals;
            vals.reserve(size);
            double_manager->getValPtr(FEM.beginKey(), FEM.endKey(), vals);
            assert(size < 7);

            int info[6];
            info[0] = size;
            for (int i = 0; i < e->size(0); ++i)
            {
               auto n = e->get(0, i);
               info[i + 1] = n->getId();
            }

            double data[16];
            for (int i = 0; i < size; ++i)
            {
               data[i] = vals[i]->getVal();
            }

            oss.write((char*)info, info_size);
            oss.write((char*)data, size * sizeof(double));
         }
         oss.close();
      }
   }
}
