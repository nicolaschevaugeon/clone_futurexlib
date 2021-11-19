/*
    This file is a part of eXlibris C++ Library
    under the GNU General Public License:
    See the LICENSE.md files for terms and
    conditions.
*/

#ifndef _xPOSTPROCESSINGMANAGER_H_
#define _xPOSTPROCESSINGMANAGER_H_

#include "xExportAlgorithm.h"
#include "xExportGmsh.h"
#include "xExportManager.h"
#include "xField.h"
#include "xSensors.h"

namespace xexport
{
class xPostProcessingManager
{
  public:
   xPostProcessingManager(xexport::xExport &, xfem::xEntityFilter = xfem::xAcceptAll());

   void initExportManager(const xParseData &, std::string, const int);
   void initSensors(xfem::xMesh &, const xParseData &, std::string, std::string);
   void reinitMesh(xfem::xMesh &);
   void setDebugFlag(bool);

   std::string getFileName() const;

   bool toExport(const std::string export_name);
   bool toExport(const std::string export_name, const int step, const std::string &extension_name = "");

   void exportOnSpace(const std::string export_name, const int step, const xfem::xLevelSet &ls);

   template <typename ITER>
   void exportOnSpace(const std::string export_name, const int step, const xfem::xEval<double> &eval,
                      const xfem::xIntegrationRule &integration_rule, ITER begin, ITER end,
                      xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity *>())
   {
      if (export_manager.toExport(export_name, step, ""))
      {
         exportOnSpace(export_name, eval, integration_rule, begin, end, integ2appro);
      }
   }

   template <typename ITER>
   void exportOnSpace(const std::string export_name, const int step, const xfem::xEval<xtensor::xVector<>> &eval,
                      const xfem::xIntegrationRule &integration_rule, ITER begin, ITER end,
                      xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity *>())
   {
      int nb_split = pexport.getNbSplit();
      if (export_manager.toExport(export_name, step, ""))
      {
         pexport.setNbSplit(nb_split);
         exportOnSpace(export_name, eval, integration_rule, begin, end, integ2appro);
      }
      for (int i = 0; i < 3; ++i)
      {
         if (export_manager.toExport(export_name + "_" + std::to_string(i + 1), step, ""))
         {
            pexport.setNbSplit(nb_split);
            xtensor::xExtractCompVector<> extract(i);
            xfem::xEvalUnary<xtensor::xExtractCompVector<>> evalex(extract, eval);
            exportOnSpace(export_name, evalex, integration_rule, begin, end, integ2appro);
         }
      }
   }

   template <typename ITER>
   void exportOnSpace(const std::string export_name, const int step, const xfem::xEval<xtensor::xTensor2<>> &eval,
                      const xfem::xIntegrationRule &integration_rule, ITER begin, ITER end,
                      xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity *>())
   {
      int nb_split = pexport.getNbSplit();
      if (export_manager.toExport(export_name, step, ""))
      {
         pexport.setNbSplit(nb_split);
         exportOnSpace(export_name, eval, integration_rule, begin, end, integ2appro);
      }
      for (int i = 0; i < 3; ++i)
      {
         for (int j = 0; j < 3; ++j)
         {
            if (export_manager.toExport(export_name + "_" + std::to_string(i + 1) + std::to_string(j + 1), step, ""))
            {
               pexport.setNbSplit(nb_split);
               xtensor::xExtractCompTensor<> extract(i, j);
               xfem::xEvalUnary<xtensor::xExtractCompTensor<>> evalex(extract, eval);
               exportOnSpace(export_name, evalex, integration_rule, begin, end, integ2appro);
            }
         }
      }
      if (export_manager.toExport(export_name + "_tr", step, ""))
      {
         pexport.setNbSplit(nb_split);
         xfem::xEvalUnary<xtensor::xTrace<>> evaltr(eval);
         exportOnSpace(export_name, evaltr, integration_rule, begin, end, integ2appro);
      }
      if (export_manager.toExport(export_name + "_dev", step, ""))
      {
         pexport.setNbSplit(nb_split);
         xfem::xEvalUnary<xtensor::xDeviatoric<>> evaldev(eval);
         exportOnSpace(export_name, evaldev, integration_rule, begin, end, integ2appro);
      }
   }

   template <typename ITER>
   void exportOnTime(const std::string name, const double time, const xfem::xEval<double> &eval,
                     const xfem::xIntegrationRule &integration_rule, ITER begin, ITER end, xfem::xEntityToEntity entity_to_entity)
   {
      if (sensor_manager.toExport(name))
      {
         double global = 0.;
         xfem::xIntegrateEvalCommand<xfem::xEval<double>> command(eval, global);
         xfem::ApplyCommandOnIntegrationRule(command, integration_rule, begin, end, entity_to_entity);
         sensor_manager.measure(global, time, name);
      }
   }

   void exportOnSpace(const std::string export_name, const int step, const xfem::xMesh &mesh);

   void exportOnTime(const std::string name, const double time, const double val);
   void exportOnTime(const std::string name, const double time, const xfem::xEval<double> &eval);
   void exportOnTime(const std::string name, const double time, const xfem::xEval<xtensor::xVector<>> &eval);
   void exportOnTime(const std::string name, const double time, const xfem::xEval<xtensor::xTensor2<>> &eval);

   void exportDebug(const std::string export_name, const int step, xfem::xValueManagerDist<double> &dm);

   void saveMesh(std::string export_name, int step, xfem::xMesh &mesh);
   void saveField(std::string export_name, int step, const xfem::xField<> &field, xfem::xRegion reg);

   template <typename ITER>
   void saveDomain(std::string export_name, int step, ITER it, ITER end)
   {
      if (step > 0 && export_manager.toExport("save", step, ""))
      {
         std::ofstream oss(export_name + "_save.bin", std::ios::binary);
         if (oss.is_open())
         {
            int nb = std::distance(it, end);
            oss.write((char *)&nb, sizeof(int));
            for (; it != end; ++it)
            {
               auto e = *it;
               int data[16];
               int size = e->size(0);
               for (int i = 0; i < size; ++i)
               {
                  data[i] = e->get(0, i)->getId();
               }
               oss.write((char *)&size, sizeof(int));
               oss.write((char *)data, size * sizeof(int));
            }
            oss.close();
         }
      }
   }

  protected:
   template <typename EVAL, typename ITER>
   void exportOnSpace(const std::string &export_name, const EVAL &eval, const xfem::xIntegrationRule &integration_rule,
                      ITER begin, ITER end, xfem::xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity *>())
   {
      if (debug)
         std::cout << "Exporting " << export_manager.getFileName() << " with nb split = " << pexport.getNbSplit() << std::endl;
      pexport.openFile(export_manager.getFileName());
      xfem::xFilteredRegion<ITER, xfem::xEntityFilter> fr(begin, end, filter);
      Export(eval, pexport, export_name, integration_rule, fr.begin(), fr.end(), integ2appro);
      pexport.closeFile();
      pexport.setNbSplit(1);
   }

   xexport::xExportManager export_manager;
   xexport::xExportSensors sensor_manager;
   xexport::xExport &pexport;
   xfem::xEntityFilter filter;
   bool debug;
};

}  // namespace xexport

#endif
