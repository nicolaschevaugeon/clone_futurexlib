/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// xfem
#include "xIntegrationRuleStored.h"

#include "xCommandOnGeomElem.h"
// xmapping
#include "xMappingBuilderHolder.h"

namespace xfem
{
using std::cout;
using std::endl;

xGaussPoints* xIntegrationRuleStoredDataManager::getStored(AOMD::mEntity& e) { return gauss_points.getData(e); }
void xIntegrationRuleStoredDataManager::delStored(AOMD::mEntity& e) { gauss_points.deleteData(e); }
xGaussPoints& xIntegrationRuleStoredDataManager::createEmptyStorageOrGet(AOMD::mEntity& e) { return gauss_points.setData(e); }
void xIntegrationRuleStoredDataManager::clearStored() { gauss_points.clear(); }

void xIntegrationRuleStoredBasic::accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const
{
   if (!filter(e_integ)) return;
   xmapping::xMapping* mapping = xMappingBuilderHolderSingleton::instance().buildMapping(*e_integ);
   xquadrature::xIntegrator* integrator = nullptr;
   const xGaussPoints* agp = gauss_points.getData(*e_integ);
   if (!agp)
   {
      integrator = new xquadrature::xGaussIntegrator(e_integ);
   }
   else
   {
      integrator = new xStoredIntegrationPoints(*agp);
   }
   xGeomElem geo_integ(e_integ, mapping, integrator);
   geo_integ.SetIntegrationPointNumberForDegree(degree);
   command.execute(&geo_integ);
   nbpoint = geo_integ.GetNbIntegrationPoints();
   delete integrator;
   delete mapping;
}

void xIntegrationRuleStoredBasicFallBackPartition::accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const
{
   const bool debug = false;

   if (!filter(e_integ)) return;
   const xGaussPoints* agp = gauss_points.getData(*e_integ);

   // If stored rule, use it
   if (agp)
   {
      if (debug) cout << "Stored\n";
      xmapping::xMapping* mapping = xMappingBuilderHolderSingleton::instance().buildMapping(*e_integ);
      xquadrature::xIntegrator* integrator = nullptr;

      integrator = new xStoredIntegrationPoints(*agp);
      xGeomElem geo_integ(e_integ, mapping, integrator);
      geo_integ.SetIntegrationPointNumberForDegree(degree);
      command.execute(&geo_integ);

      nbpoint = geo_integ.GetNbIntegrationPoints();
      delete integrator;
      delete mapping;
   }
   else
   {
      // Else use partition
      if (debug) cout << "FallBack\n";
      xPartition partition;
      getpartition(e_integ, partition, filter);

      for (AOMD::mEntity* es_integ : partition)
      {
         if (debug)
         {
            std::cout << "partition found inside sub\n";
            es_integ->print();
         }
         xmapping::xMapping* mapping = xMappingBuilderHolderSingleton::instance().buildMapping(*es_integ);
         xGeomElem geo_integ_es(es_integ, mapping, getIntegrator(es_integ));
         geo_integ_es.SetIntegrationPointNumberForDegree(degree);
         nbpoint += geo_integ_es.GetNbIntegrationPoints();
         command.execute(&geo_integ_es);
         delete mapping;
      }
   }
   return;
}

void xIntegrationRuleStoredPartition::accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const
{
   const bool debug = false;

   xPartition partition;
   if (debug) std::cout << "On ouvre l element ";
   if (debug) e_integ->print();
   getpartition(e_integ, partition, filter);

   if (debug) std::cout << "partition size is " << partition.size() << std::endl;

   for (AOMD::mEntity* es_integ : partition)
   {
      if (debug) std::cout << "partition found inside sub\n";
      if (debug) es_integ->print();
      xmapping::xMapping* mapping = xMappingBuilderHolderSingleton::instance().buildMapping(*es_integ);
      xquadrature::xIntegrator* integrator = nullptr;
      const xGaussPoints* agp = gauss_points.getData(*es_integ);
      if (!agp)
      {
         integrator = new xquadrature::xGaussIntegrator(es_integ);
      }
      else
      {
         integrator = new xStoredIntegrationPoints(*agp);
      }
      xGeomElem geo_integ_es(es_integ, mapping, integrator);
      geo_integ_es.SetIntegrationPointNumberForDegree(degree);
      nbpoint += geo_integ_es.GetNbIntegrationPoints();
      command.execute(&geo_integ_es);
      delete integrator;
      delete mapping;
   }
   if (debug) std::cout << "--------\n";
}

}  // namespace xfem
