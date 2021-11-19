/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <cmath>
#include <iostream>
#include <vector>
// aomd
#include "mEntity.h"

// xfem
#include "xCommandOnGeomElem.h"
#include "xElement.h"
#include "xEntityToEntity.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xLevelSet.h"
#include "xMappingBuilderHolder.h"

using AOMD::mEntity;

namespace xfem
{
std::vector<xtensor::xPoint>* rectangleRulePointsTri[xRectangleIntegrator::levelmax] = {nullptr, nullptr, nullptr, nullptr,
                                                                                        nullptr, nullptr, nullptr, nullptr};

int xRectangleIntegrator::nbIntegrationPoints(int level) const { return static_cast<int>(std::pow(4., level)); }

void xRectangleIntegrator::iPoint(int i, int level, double& u, double& v, double& w, double& weight) const
{
   if (level > levelmax) throw;
   // only triangle for now ... add a switch
   if (!rectangleRulePointsTri[level])
   {
      rectangleRulePointsTri[level] = new std::vector<xtensor::xPoint>(integpoint(level));
   }
   xtensor::xPoint& P = (*rectangleRulePointsTri[level])[i];
   u = P(0);
   v = P(1);
   w = P(2);
   weight = 0.5 / rectangleRulePointsTri[level]->size();
}

xIntegrationRule::xIntegrationRule() : currentRule(GAUSS), integrator(nullptr) {}

void xIntegrationRuleBasic::accept(xCommandOnGeomElem& command, mEntity* e_integ) const
{
   if (!filter(e_integ)) return;
   xGeomElem geo_integ(e_integ);
   geo_integ.SetIntegrationPointNumberForDegree(degree);
   command.execute(&geo_integ);
   nbpoint = geo_integ.GetNbIntegrationPoints();
}

void xIntegrationRuleNodal::accept(xCommandOnGeomElem& command, mEntity* e_integ) const
{
   const bool debug = false;
   if (!filter(e_integ)) return;
   xmapping::xMapping* mapping = xMappingBuilderHolderSingleton::instance().buildMapping(*e_integ);
   xNodalIntegrationPoints integration_points(e_integ);
   xGeomElem geo_integ(e_integ, mapping, &integration_points);
   geo_integ.SetIntegrationPointNumberForDegree(0);
   if (debug) std::cout << "nb of integration points " << geo_integ.GetNbIntegrationPoints() << std::endl;
   command.execute(&geo_integ);
   delete mapping;
}

void xIntegrationRulePartition::accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const
{
   // nbpoint = 0;
   const bool debug = false;

   xPartition partition;
   if (debug) std::cout << "On ouvre l element ";
   if (debug) e_integ->print();
   getpartition(e_integ, partition, filter);

   if (debug) std::cout << "partition size is " << partition.size() << std::endl;

   for (mEntity* es_integ : partition)
   {
      if (debug) std::cout << "partition found inside sub\n";
      if (debug) es_integ->print();
      xmapping::xMapping* mapping = xMappingBuilderHolderSingleton::instance().buildMapping(*es_integ);
      xGeomElem geo_integ_es(es_integ, mapping, getIntegrator(es_integ));
      geo_integ_es.SetIntegrationPointNumberForDegree(degree);
      nbpoint += geo_integ_es.GetNbIntegrationPoints();
      command.execute(&geo_integ_es);
      delete mapping;
   }
   if (debug) std::cout << "--------\n";
}

void xIntegrationRulePartitionBoundary::accept(xCommandOnGeomElem& command, mEntity* e_integ) const
{
   const bool debug = false;
   if (debug) std::cout << "On ouvre l element ";
   if (debug) e_integ->print();

   xUpperAdjacencyPartition upperget;
   xPartition partition;
   getpartition(e_integ, partition, xAcceptAll());
   if (debug) cout << "partition size is " << partition.size() << std::endl;
   for (mEntity* es_integ : partition)
   {
      if (debug)
      {
         std::cout << "partition found inside sub for xIntegrationRulePartitionBoundary \n";
         std::cout << "On ouvre l element es_integ";
         es_integ->print();
      }

      if ((e_integ == es_integ) && (filter(upperget(es_integ))))
      {
         if (debug) cout << " test ((e_integ == es_integ)&&(filter(upperget(es_integ)))) VRAI" << endl;
         xGeomElem geo_integ(e_integ);
         geo_integ.SetIntegrationPointNumberForDegree(degree);
         command.execute(&geo_integ);
      }

      else if ((e_integ != es_integ) && (filter(es_integ)))
      {
         if (debug) cout << " test ((e_integ != es_integ)&&(filter(es_integ))) VRAI " << endl;
         xGeomElem geo_integ_es(es_integ);
         geo_integ_es.SetIntegrationPointNumberForDegree(degree);
         command.execute(&geo_integ_es);
      }
   }
}

void xIntegrationRulePartitionExport::accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const
{
   const bool debug = false;
   xPartition partition;
   if (debug) std::cout << "On ouvre l element ";
   if (debug) e_integ->print();

   std::vector<double> c = lst->getVals(e_integ);

   // moche : a revoir...
   double min = c[0], max = c[0];
   for (size_t i = 1; i < c.size(); ++i)
   {
      if (c[i] > max) max = c[i];
      if (c[i] < min) min = c[i];
   }

   if (debug) cout << "min=" << min << " max=" << max << endl;

   if (min * max < 0. || max < 0.)
   {
      getpartition(e_integ, partition, filter);
      if (debug) std::cout << "partition size is " << partition.size() << std::endl;
      for (mEntity* es_integ : partition)
      {
         if (debug)
         {
            std::cout << "partition found inside sub\n";
            es_integ->print();
         }
         xGeomElem geo_integ_es(es_integ);
         geo_integ_es.SetIntegrationPointNumberForDegree(degree);
         command.execute(&geo_integ_es);
      }
      if (debug) std::cout << "--------\n";
   }
   else
   {
      if (!filter(e_integ)) return;
      xGeomElem geo_integ(e_integ);
      geo_integ.SetIntegrationPointNumberForDegree(degree);
      command.execute(&geo_integ);
   }
}

xNodalIntegrationPoints::xNodalIntegrationPoints(pEntity e)
//  : ent(e)
{
   // entityType = M_GetElementType(e);
   entityType = e->getType();
}

int xNodalIntegrationPoints::nbIntegrationPoints(int order) const
{
   // correction made by N. MOES at Ecole Centrale Nantes
   switch (entityType)
   {
      case AOMD::mEntity::EDGE:
      {
         return 2;
         break;
      }
      case AOMD::mEntity::TRI:
      {
         return 3;
         break;
      }
      default:
         std::cerr << "not coded yet" << std::endl;
         assert(0);
         throw;
         break;
   }
}
void xNodalIntegrationPoints::iPoint(int i, int order, double& u, double& v, double& w, double& weight) const
{
   switch (entityType)
   {
      case AOMD::mEntity::EDGE:
      {
         {
            if (i == 0)
               u = -1.0;
            else if (i == 1)
               u = 1.0;
            else
               assert(0);
            v = w = 0.0;
            weight = 1.;
         }
         break;
      }
      case AOMD::mEntity::TRI:
      {
         switch (i)
         {
            case 0:
            {
               u = 0.;
               v = 0.;
               w = 0.;
               weight = 1. / 6.;
               break;
            }
            case 1:
            {
               u = 1.;
               v = 0.;
               w = 0.;
               weight = 1. / 6.;
               break;
            }
            case 2:
            {
               u = 0.;
               v = 1.;
               w = 0.;
               weight = 1. / 6.;
               break;
            }
         }
         break;
      }
      default:
         std::cerr << "not coded yet" << std::endl;
         assert(0);
         break;
   }
}

/// K : for integration of discontinuous functions

void integpointnewlevelG(int u0, int v0, int* u, int* v)
{
   int u0n = u0 * 2;
   int v0n = v0 * 2;

   u[0] = u0n - 1;
   v[0] = v0n - 1;

   u[1] = u0n + 2;
   v[1] = v0n - 1;

   u[2] = u0n - 1;
   v[2] = v0n + 2;

   u[3] = u0n;
   v[3] = v0n;
}

void integpointnewlevelD(int u0, int v0, int* u, int* v)
{
   int u0n = u0 * 2;
   int v0n = v0 * 2;

   u[0] = u0n + 1;
   v[0] = v0n - 2;

   u[1] = u0n + 1;
   v[1] = v0n + 1;

   u[2] = u0n - 2;
   v[2] = v0n + 1;

   u[3] = u0n;
   v[3] = v0n;
}

//

std::vector<xtensor::xPoint> integpoint(int level)
{
   int nptf = static_cast<int>(std::pow(4., level));
   int denom = 3 * static_cast<int>(std::pow(2., level));

   std::vector<int> uOld(1), vOld(1);
   std::vector<int> uNew(1), vNew(1);
   std::vector<bool> DGOld(1), DGNew(1);
   uNew[0] = 1, vNew[0] = 1;
   DGNew[0] = 1;

   //  res[0] = xtensor::xPoint (((double)u0)/k, ((double)v0)/k , 0.);

   for (int i = 0; i < level; ++i)
   {
      int npt = static_cast<int>(std::pow(4., i + 1));
      uOld = uNew;
      vOld = vNew;
      DGOld = DGNew;
      uNew.clear();
      uNew.resize(npt);
      vNew.clear();
      vNew.resize(npt);
      DGNew.clear();
      DGNew.resize(npt);
      for (unsigned int k = 0; k < uOld.size(); ++k)
      {
         if (DGOld[k])
         {
            integpointnewlevelG(uOld[k], vOld[k], &uNew[k * 4], &vNew[k * 4]);
            DGNew[k * 4] = 1;
            DGNew[k * 4 + 1] = 1;
            DGNew[k * 4 + 2] = 1;
            DGNew[k * 4 + 3] = 0;
         }
         else
         {
            integpointnewlevelD(uOld[k], vOld[k], &uNew[k * 4], &vNew[k * 4]);
            DGNew[k * 4] = 0;
            DGNew[k * 4 + 1] = 0;
            DGNew[k * 4 + 2] = 0;
            DGNew[k * 4 + 3] = 1;
         }
      }
   }

   std::vector<xtensor::xPoint> res(nptf);
   for (int i = 0; i < nptf; ++i)
   {
      res[i] = xtensor::xPoint(((double)uNew[i]) / denom, ((double)vNew[i]) / denom, 0.);
   }
   return res;
}

}  // namespace xfem
