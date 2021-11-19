#include "xApproxFunctionBySupportComponent.h"

#include "mVertex.h"
#include "xElement.h"
#include "xGeomElem.h"
#include "xMesh.h"

namespace xfem
{
xHeavisideBySub::xHeavisideBySub(std::unordered_set<const AOMD::mEntity *> _data) : data(std::move(_data)) {}

void xHeavisideBySub::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   const AOMD::mEntity *e = geo_integ->getEntity();
   res = (data.count(e)) ? 1. : 0.;
}

void xHeavisideBySub::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ, xtensor::xVector<> &res) const
{
   res = {0., 0., 0.};
}

xLinearBySub::xLinearBySub(std::unordered_map<const AOMD::mVertex *, double> _data) : data(std::move(_data)) {}
void xLinearBySub::getVal(const xGeomElem *geo_appro, const xGeomElem *geo_integ, double &res) const
{
   const AOMD::mEntity *e = geo_integ->getEntity();
   xfem::xElement elem(e);
   elem.setUvw(geo_integ->getUVW());
   res = elem.getInterpoSca(vals(*e));
}
void xLinearBySub::getGrad(const xGeomElem *geo_appro, const xGeomElem *geo_integ, xtensor::xVector<> &res) const
{
   const AOMD::mEntity *e = geo_integ->getEntity();
   xfem::xElement elem(e);
   elem.setUvw(geo_integ->getUVW());
   res = elem.getGradInterpoSca(vals(*e));
}

std::vector<double> xLinearBySub::vals(const AOMD::mEntity &e_integ) const
{
   const size_t nb_vertex = size_t(e_integ.size(0));
   std::vector<double> vals(nb_vertex);
   for (size_t i = 0; i < nb_vertex; ++i)
   {
      const AOMD::mVertex *v = static_cast<const AOMD::mVertex *>(e_integ.get(0, int(i)));
      auto it = data.find(v);
      if (it != data.end())
         vals[i] = it->second;
      else
         vals[i] = 0.;
   }
   return vals;
}

xApproxFunctionBySupportComponentGenerator::xApproxFunctionBySupportComponentGenerator(
    std::function<const xfem::xSupportComponent *(AOMD::mEntity *e)> _getSupportComponent)
    : getSupportComponent(_getSupportComponent)
{
}

void xApproxFunctionBySupportComponentGenerator::generateKeysAndFcts(AOMD::mEntity *e, xfem::xValKey &key,
                                                                     xfem::ValKeyModifierPtr key_modifier, approxFctPtr_t f,
                                                                     femKeys *keys, femFcts_t *fcts)
{
   AOMD::mEntity *e_enriched = key.getEnti();
   femFcts_t fctsnew;
   getEnrichmentFunctions(e_enriched, fctsnew);
   for (size_t i = 0; i < fctsnew.size(); ++i)
   {
      fcts->push_back(approxFctPtr_t(new xApproxFunctionEnrichedXFEM(f, fctsnew[i])));
      key_modifier->addToExtension(std::to_string(i));
      (*key_modifier)(key);
      keys->push_back(key);
   }
}
void xApproxFunctionBySupportComponentGenerator::generateKeys(xfem::xValKey &key, xfem::ValKeyModifierPtr key_modifier,
                                                              femKeys *keys)
{
   AOMD::mEntity *e_enriched = key.getEnti();
   femFcts_t fctsnew;
   getEnrichmentFunctions(e_enriched, fctsnew);
   for (size_t i = 0; i < fctsnew.size(); ++i)
   {
      key_modifier->addToExtension(std::to_string(i));
      (*key_modifier)(key);
      keys->push_back(key);
   }
}
void xApproxFunctionBySupportComponentGenerator::clear() { generated_functions.clear(); }
xPartition xApproxFunctionBySupportComponentGenerator::getPartition(const AOMD::mEntity &e) const
{
   xPartition part;
   xMesh::getPartition(const_cast<AOMD::mEntity *>(&e), part, xAcceptAll());
   return part;
}

void xApproxFunctionConstantBySupportComponentGenerator::getEnrichmentFunctions(AOMD::mEntity *enrichedentity,
                                                                                femFcts_t &functions)
{
   functions.clear();
   femFcts_t *pfunctions = generated_functions.getData(*enrichedentity);
   if (pfunctions)
   {
      functions = *pfunctions;
      return;
   }
   if (!pfunctions)
   {
      const xfem::xSupportComponent *psupportcomponent = getSupportComponent(enrichedentity);
      if (!psupportcomponent) return;
      const xfem::xSupportComponent &supportcomponent = *psupportcomponent;
      const size_t nbcomponent = supportcomponent.getnbComponents();
      if (nbcomponent < 2) return;
      functions.resize(nbcomponent - 1);
      for (size_t i = 0; i < nbcomponent - 1; ++i)
      {
         const xfem::xSupportComponent::component &componenti = supportcomponent(i);
         std::unordered_set<const AOMD::mEntity *> entities_at_one;
         for (const AOMD::mEntity *e : componenti)
         {
            entities_at_one.insert(e);
            // e may be itself cut by a level set : => tag partition
            for (const AOMD::mEntity *e_part : getPartition(*e)) entities_at_one.insert(e_part);
         }
         functions[i] = approxFctPtr_t(new xHeavisideBySub(entities_at_one));
      }
      generated_functions.setData(*enrichedentity) = functions;
   }
}

void xApproxFunctionConstantAndRampBySupportComponentGenerator::getEnrichmentFunctions(AOMD::mEntity *enrichedentity,
                                                                                       femFcts_t &functions)
{
   functions.clear();
   femFcts_t *pfunctions = generated_functions.getData(*enrichedentity);
   if (pfunctions)
   {
      functions = *pfunctions;
      return;
   }
   if (!pfunctions)
   {
      const xfem::xSupportComponent *psupportcomponent = getSupportComponent(enrichedentity);
      if (!psupportcomponent) return;
      const xfem::xSupportComponent &supportcomponent = *psupportcomponent;
      const size_t nbcomponent = supportcomponent.getnbComponents();
      if (nbcomponent < 2) return;
      functions.resize(nbcomponent - 1);
      for (size_t i = 0; i < nbcomponent - 1; ++i)
      {
         std::unordered_map<const AOMD::mVertex *, double> vertex_values;
         const xfem::xSupportComponent::component &ramp_part = supportcomponent.getSplitZone();
         for (const AOMD::mEntity *e : ramp_part)
         {
            for (const AOMD::mEntity *e_part : getPartition(*e))
            {
               for (size_t iv = 0; iv < size_t(e_part->size(0)); ++iv)
               {
                  const AOMD::mVertex *v = static_cast<const AOMD::mVertex *>(e_part->get(0, int(iv)));
                  vertex_values[v] = 1. / (nbcomponent);
               }
            }
         }
         for (size_t k = 0; k < nbcomponent; ++k)
         {
            if (k != i)
            {
               const xfem::xSupportComponent::component &componentk = supportcomponent(k);
               for (const AOMD::mEntity *e : componentk)
               {
                  for (const AOMD::mEntity *e_part : getPartition(*e))
                  {
                     for (size_t iv = 0; iv < size_t(e_part->size(0)); ++iv)
                     {
                        const AOMD::mVertex *v = static_cast<const AOMD::mVertex *>(e_part->get(0, int(iv)));
                        vertex_values[v] = 0.;
                     }
                  }
               }
            }
         }
         for (const AOMD::mEntity *e : supportcomponent(i))
         {
            for (const AOMD::mEntity *e_part : getPartition(*e))
            {
               for (size_t iv = 0; iv < size_t(e_part->size(0)); ++iv)
               {
                  const AOMD::mVertex *v = static_cast<const AOMD::mVertex *>(e_part->get(0, int(iv)));
                  vertex_values[v] = 1.;
               }
            }
         }
         functions[i] = approxFctPtr_t(new xLinearBySub(vertex_values));
      }
      generated_functions.setData(*enrichedentity) = functions;
   }
}

}  // namespace xfem
