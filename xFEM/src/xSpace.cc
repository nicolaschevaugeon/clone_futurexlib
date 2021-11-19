/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xSpace.h"

#include <string>
#include <vector>

#include "xApproxFunction.h"
#include "xApproxFunctionCorrectedShifted.h"
#include "xApproxFunctionHighOrder.h"
#include "xFiniteElement.h"
#include "xIntegrationRule.h"
#include "xValKey.h"

using AOMD::mEntity;
using std::cout;
using std::endl;
namespace xfem
{
xSpaceRegular::xSpaceRegular(const std::string& a, TensorialType_t b) : Phys(xKeyInfo::getPhysId(a)), TensorialType(b) {}

xSpacePointwise::xSpacePointwise(xIntegrationRule& integ_) : integ(integ_), GP_KEY(2000)  // max 2000 Gauss points per element !
{
   for (unsigned int i = 0; i < GP_KEY.size(); ++i)
   {
      char str[10];
      std::sprintf(str, "%u", i + 1);
      std::string name("GAUSS_" + std::string(str));
      GP_KEY[i] = xKeyInfo::getGeomId(name);
   }
}

void xSpacePointwise::set_GP_KEY_size(unsigned int limit)
{
   GP_KEY.resize(limit, 0);
   for (unsigned int i = 0; i < GP_KEY.size(); ++i)
   {
      char str[10];
      std::sprintf(str, "%u", i + 1);
      std::string name("GAUSS_" + std::string(str));
      GP_KEY[i] = xKeyInfo::getGeomId(name);
   }
}

void xSpacePwRegular::getKeys(mEntity* e, femKeys* keys)
{
   xGetKeysOnGeomElemCommand command(keys, GP_KEY, Phys);
   xGeomElem geom_appro(e);
   command.openApproxElem(&geom_appro);
   integ.accept(command, e);
   command.closeApproxElem(&geom_appro);
}

void xSpacePwRegular::getKey(const xGeomElem* geo_integ, xValKey& key)
{
   key = xValKey(Phys, GP_KEY[geo_integ->getCurrentIntegrationPointId()], geo_integ->getEntity());
}

xSpaceFiltered::~xSpaceFiltered() = default;

void xSpaceFiltered::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   femKeys keysl;
   femFcts fctsl;
   base_space->getKeysAndFcts(e, &keysl, &fctsl);
   // filtrage
   auto itfl = fctsl.begin();
   for (auto itkl = keysl.begin(); itkl != keysl.end(); ++itkl, ++itfl)
   {
      if (accept(itkl->getEnti()))
      {
         keys->push_back(*itkl);
         fcts->push_back(*itfl);
      }
   }
   return;
   /*
   //  old implementation : erase the key and function : pb filter is applied on all the keys, not the one just got from the base
   space of the filter space. base_space->getKeysAndFcts(e, keys, fcts); femKeys::iterator  it = keys->begin(); femFcts::iterator
   itf = fcts->begin(); for (; it != keys->end(); )
     {
       mEntity* e = it->getEnti();
       if (!accept(e))
         {
           it  = keys->erase(it);
           itf = fcts->erase(itf);
         }
       else
         {
           ++it;
           ++itf;
         }
     }
   return;
   */
}

void xSpaceFiltered::getKeys(mEntity* e, femKeys* keys)
{
   femKeys keysl;
   base_space->getKeys(e, &keysl);
   // filtrage
   std::copy_if(keysl.begin(), keysl.end(), std::back_inserter(*keys),
                [this](const xValKey& key) { return accept(key.getEnti()); });
   return;

   /*
   //  old implementation : erase the key and function : pb filter is applied on all the keys, not the one just got from the base
   space of the filter space. base_space->getKeys(e, keys);
    //filtrage
    femKeys::iterator it = keys->begin();
    for (; it != keys->end(); )
      {
        mEntity* e = it->getEnti();
        if (!accept(e)) it = keys->erase(it);
        else ++it;
      }

   return;
   */
}

xSpaceKeyFiltered::~xSpaceKeyFiltered() = default;

void xSpaceKeyFiltered::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   femKeys keysl;
   femFcts fctsl;
   base_space->getKeysAndFcts(e, &keysl, &fctsl);
   auto itfl = fctsl.begin();
   for (auto itkl = keysl.begin(); itkl != keysl.end(); ++itkl, ++itfl)
   {
      if (accept(*itkl))
      {
         keys->push_back(*itkl);
         fcts->push_back(*itfl);
      }
   }
   return;
   // filtrage

   // old implementation : wrong : suppose the input keys and fcts are empty.
   /*
   base_space->getKeysAndFcts(e, keysl, fctsl);
   femKeys::iterator  it = keys->begin();
   femFcts::iterator itf = fcts->begin();
         for (; it != keys->end(); )
         {
                 if (!accept((*it)))
                 {
                         it  = keys->erase(it);
                         itf = fcts->erase(itf);
                 }
                 else
                 {
                         ++it;
                         ++itf;
                 }
         }
         return;
   */
}

void xSpaceKeyFiltered::getKeys(mEntity* e, femKeys* keys)
{
   femKeys keysb;
   base_space->getKeys(e, &keysb);
   std::copy_if(keysb.begin(), keysb.end(), std::back_inserter(*keys), [this](xValKey& key) { return accept(key); });
   return;
   // old implementation buggy : suppose the keys arrive empty.
   /*
         femKeys::iterator it = keys->begin();
         for (; it != keys->end(); )
         {
                 if (!accept((*it))) it = keys->erase(it);
                 else ++it;
         }

         return;
   */
}

xDiscontinuousSpaceFiltered::~xDiscontinuousSpaceFiltered() = default;

void xDiscontinuousSpaceFiltered::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   if (accept(e)) base_space->getKeysAndFcts(e, keys, fcts);
   return;
}

void xDiscontinuousSpaceFiltered::getKeys(mEntity* e, femKeys* keys)
{
   if (accept(e)) base_space->getKeys(e, keys);
   return;
}

void xSpaceComposite::getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* fcts)
{
   for (const_iterator it = children.begin(); it != children.end(); ++it)
   { /*
      femKeys loc;
      femFcts fct;
      (*it)->getKeysAndFcts(e, &loc, &fct);
      keys->insert(keys->end(), loc.begin(), loc.end());
      fcts->insert(fcts->end(), fct.begin(), fct.end());
     */
      (*it)->getKeysAndFcts(e, keys, fcts);
   }

   return;
}

void xSpaceComposite::getKeys(AOMD::mEntity* e, femKeys* keys)
{
   for (const_iterator it = children.begin(); it != children.end(); ++it)
   {
      /*
      femKeys loc;
      (*it)->getKeys(e, &loc);
      keys->insert(keys->end(), loc.begin(), loc.end());
      */
      (*it)->getKeys(e, keys);
   }

   return;
}

void xSpaceDifference::getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* fcts)
{
   femKeys first_vector_key, second_vector_key;
   femFcts first_vector_fct, second_vector_fct;
   first->getKeysAndFcts(e, &first_vector_key, &first_vector_fct);
   second->getKeysAndFcts(e, &second_vector_key, &second_vector_fct);

   // implementation not efficient if first and second are big
   std::vector<xValKey>::iterator it;
   for (it = second_vector_key.begin(); it != second_vector_key.end(); ++it)
   {
      std::vector<xValKey>::iterator itk = find(first_vector_key.begin(), first_vector_key.end(), *it);
      if (itk != first_vector_key.end())
      {
         femFcts::iterator itf = first_vector_fct.begin();
         advance(itf, distance(first_vector_key.begin(), itk));
         first_vector_key.erase(itk);
         first_vector_fct.erase(itf);
      }
   }
   keys->insert(keys->end(), first_vector_key.begin(), first_vector_key.end());
   fcts->insert(fcts->end(), first_vector_fct.begin(), first_vector_fct.end());
   return;
}

void xSpaceDifference::getKeys(AOMD::mEntity* e, femKeys* keys)
{
   femKeys first_vector, second_vector;
   first->getKeys(e, &first_vector);
   second->getKeys(e, &second_vector);

   // implementation not efficient if first and second are big
   std::vector<xValKey>::iterator it;
   for (it = second_vector.begin(); it != second_vector.end(); ++it)
   {
      first_vector.erase(std::remove(first_vector.begin(), first_vector.end(), *it), first_vector.end());
   }
   keys->insert(keys->end(), first_vector.begin(), first_vector.end());

// set implementation to improve speed when first and second are big
// not debugged yet
#if 0
  std::set<xValKey> first_set(first_vector.begin(), first_vector.second());
  std::set<xValKey> second_set(first_vector.begin(), first_vector.second());
  iterKey end = std::set_difference(first_set.begin(),  first_set.end(),
				    second_set.begin(), second_set.end(),
				    first_vector.begin());
  first_vector.erase(end,  first_vector.end());
  keys->insert(keys->end(), first_vector.begin(), fist_vector.end());
#endif

   return;
}

void xSpaceConstant::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   int first = keys->size();
   getKeys(e, keys);
   int last = keys->size();
   switch (TensorialType)
   {
      case SCALAR:
         for (int i = first; i < last; ++i) appro->push_back(shapeFctPtr(new xApproxFunctionConstant()));
         break;
      default:
         assert(0);
         break;
   }
}

void xSpaceConstant::getKeys(mEntity* e, femKeys* keys)
{
   keys->push_back(xValKey(Phys, CONSTANT_ELEMENT, e));
   return;
}

void xSpaceConstantGlobal::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   int first = keys->size();
   getKeys(global_entity, keys);
   int last = keys->size();
   switch (TensorialType)
   {
      case SCALAR:
         for (int i = first; i < last; ++i) appro->push_back(shapeFctPtr(new xApproxFunctionConstant()));
         break;
      default:
         assert(0);
         break;
   }
}

void xSpaceConstantGlobal::getKeys(mEntity* e, femKeys* keys)
{
   keys->push_back(xValKey(Phys, CONSTANT_ELEMENT, global_entity));
   return;
}

void xSpaceLagrange::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   // warnign if e is a vertex this function does not have to much sense
   size_t start = keys->size();
   getKeys(e, keys);
   femKeys::iterator it = keys->begin() + start;
   for (int j = 1; it != keys->end(); ++it, ++j)
   {
      switch (TensorialType)
      {
         case SCALAR:
            appro->push_back(shapeFctPtr(new xApproxFunctionScalarHierarchical(j)));
            break;
         case VECTOR_X:
            appro->push_back(shapeFctPtr(new xApproxFunctionVectorHierarchical(j, 0)));
            break;
         case VECTOR_Y:
            appro->push_back(shapeFctPtr(new xApproxFunctionVectorHierarchical(j, 1)));
            break;
         case VECTOR_Z:
            appro->push_back(shapeFctPtr(new xApproxFunctionVectorHierarchical(j, 2)));
            break;
      }
   }
}

void xSpaceLagrange::getKeys(mEntity* e, femKeys* keys)
{
   int i;
   // std::cout << "in getKeys. level : " <<  e->getLevel() << std::endl;
   if (e->getLevel() == 0)
   {
      keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e));
   }
   else if (e->getLevel() == 1)
   {
      switch (Degree)
      {
         case DEGREE_ONE:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            break;
         case DEGREE_TWO:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e));
            break;
         case DEGREE_THREE:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e));
            keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[1], e));
            break;
         default:
            assert(1 == 0);
            break;
      }
   }
   else if (e->getLevel() == 2)
   {
      switch (Degree)
      {
         case DEGREE_ONE:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            break;
         case DEGREE_TWO:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            for (i = 0; i < e->size(1); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(1, i)));
            }
            if (e->getType() == mEntity::QUAD)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e));
            }
            break;
         case DEGREE_THREE:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            for (i = 0; i < e->size(1); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(1, i)));
            }
            for (i = 0; i < e->size(1); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[1], e->get(1, i)));
            }
            if (e->getType() == mEntity::TRI)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e));
            }
            if (e->getType() == mEntity::QUAD)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e));
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[1], e));
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[2], e));
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[3], e));
            }
            break;
         default:
            assert(1 == 0);
            break;
      }
   }
   //
   // Dimension 3
   //
   else if (e->getLevel() == 3)
   {
      switch (Degree)
      {
         case DEGREE_ONE:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            break;
         case DEGREE_TWO:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            for (i = 0; i < e->size(1); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(1, i)));
            }
            break;
         case DEGREE_THREE:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            for (i = 0; i < e->size(1); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(1, i)));
            }
            for (i = 0; i < e->size(1); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[1], e->get(1, i)));
            }

            break;
         default:
            assert(1 == 0);
            break;
      }
   }
   else
   {
      fprintf(stderr, "wrong Dimension \n");
      assert(1 == 0);
   }

   return;
}
void xSpacePieceWiseLagrange::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* appro)
{
   // warnign if e is a vertex this function does not have to much sense
   size_t start = keys->size();
   getKeys(e, keys);
   femKeys::iterator it = keys->begin() + start;
   for (int j = 1; it != keys->end(); ++it, ++j)
   {
      switch (TensorialType)
      {
         case SCALAR:
            appro->push_back(shapeFctPtr(new xApproxFunctionScalarPieceWiseHierarchical(j)));
            break;
         default:
            assert(0);  // to be coded;
            break;
      }
   }
}

void xSpacePieceWiseLagrange::getKeys(mEntity* e, femKeys* keys)
{
   int i;
   // std::cout << "in getKeys. level : " <<  e->getLevel() << std::endl;
   if (e->getLevel() == 0)
   {
      keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e));
   }
   else if (e->getLevel() == 1)
   {
      switch (Degree)
      {
         case DEGREE_ZERO:
            for (i = 0; i < e->size(0); ++i)
            {
               keys->push_back(xValKey(Phys, HIERARCHICAL_KEY[0], e->get(0, i)));
            }
            break;
         default:
            assert(1 == 0);
            break;  // to be coded;
      }
   }
   else if (e->getLevel() == 2)
   {
      assert(1 == 0);  // to be coded
   }
   else if (e->getLevel() == 3)
   {
      assert(1 == 0);  // to be coded
   }
   else
   {
      fprintf(stderr, "wrong Dimension \n");
      assert(1 == 0);
   }

   return;
}

void xSpaceXFEM::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   size_t start = keys->size();
   base_space->getKeysAndFcts(e, keys, fcts);
   femKeys::iterator it = keys->begin() + start;
   femFcts::iterator itf = fcts->begin() + start;
   for (; it != keys->end(); ++it, ++itf)
   {
      *itf = shapeFctPtr(new xApproxFunctionEnrichedXFEM(*itf, enrichment));
      key_modifier(*it);
   }
   return;
}

void xSpaceVectorXFEM::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   size_t start = keys->size();
   base_space->getKeysAndFcts(e, keys, fcts);
   femKeys::iterator it = keys->begin() + start;
   femFcts::iterator itf = fcts->begin() + start;
   for (; it != keys->end(); ++it, ++itf)
   {
      *itf = shapeFctPtr(new xApproxFunctionEnrichedVectorXFEM(*itf, enrichment));
      key_modifier(*it);
   }
   return;
}

void xSpaceVector2XFEM::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   size_t start = keys->size();
   base_space->getKeysAndFcts(e, keys, fcts);
   femKeys::iterator it = keys->begin() + start;
   femFcts::iterator itf = fcts->begin() + start;
   for (; it != keys->end(); ++it, ++itf)
   {
      // *itf = shapeFctPtr(new xApproxFunctionEnrichedVector2XFEM(*itf, enrichment, it->getEnti()));
      *itf = shapeFctPtr(new xApproxFunctionEnrichedVector2XFEM(*itf, enrichment, *it));

      key_modifier(*it);
   }
   return;
}

void xSpaceXFEM::getKeys(mEntity* e, femKeys* keys)
{
   size_t start = keys->size();
   base_space->getKeys(e, keys);
   for_each(keys->begin() + start, keys->end(), key_modifier);
   return;
}

void xSpaceXFEMShifted::getKeys(mEntity* e, femKeys* keys)
{
   size_t start = keys->size();
   base_space->getKeys(e, keys);
   for_each(keys->begin() + start, keys->end(), key_modifier);
   return;
}

void xSpaceXFEMShifted::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   size_t start = keys->size();
   base_space->getKeysAndFcts(e, keys, fcts);
   femKeys::iterator it = keys->begin() + start;
   femFcts::iterator itf = fcts->begin() + start;
   for (; it != keys->end(); ++it, ++itf)
   {
      *itf = shapeFctPtr(new xApproxFunctionEnrichedXFEMShifted(*itf, enrichment, it->getEnti()));
      key_modifier(*it);
   }
   return;
}

void xSpaceXFEMShiftedCorrected::getKeys(mEntity* e, femKeys* keys)
{
   size_t start = keys->size();
   base_space->getKeys(e, keys);
   for_each(keys->begin() + start, keys->end(), key_modifier);
   return;
}

void xSpaceXFEMShiftedCorrected::getKeysAndFcts(mEntity* e, femKeys* keys, femFcts* fcts)
{
   size_t start = keys->size();
   base_space->getKeysAndFcts(e, keys, fcts);
   femKeys::iterator it = keys->begin() + start;
   femFcts::iterator itf = fcts->begin() + start;

   for (; it != keys->end(); ++it, ++itf)
   {
      shapeFctPtr scalfct =
          dynamic_cast<xApproxFunctionVector*>(dynamic_cast<xApproxFunctionVector*>((*itf).get()))->getScalarFunc();
      xtensor::xPoint uvw = getInterpolantUVW(*it, scalfct);
      *itf = shapeFctPtr(new xApproxFunctionEnrichedXFEMShiftedCorr(*itf, enrichment, it->getEnti(), uvw, ramp, gramp));
      key_modifier(*it);
   }
   return;
}

xtensor::xPoint xSpaceXFEMShiftedCorrected::getInterpolantUVW(xValKey& key, shapeFctPtr& func)
{
   mEntity* enti = key.getEnti();

   // Point
   if (enti->getLevel() == 0)
   {
      return xtensor::xPoint(0., 0., 0.);
   }

   // Edge
   if (enti->getLevel() == 1)
   {
      xApproxFunctionScalarEdgeLagrange* shapeFuncE = dynamic_cast<xApproxFunctionScalarEdgeLagrange*>((func.get()));
      if (shapeFuncE)
      {
         std::vector<double> nodeuvw;
         shapeFuncE->getNodeUvw(nodeuvw);
         return xtensor::xPoint(nodeuvw[0], nodeuvw[1], nodeuvw[2]);
      }
      else
      {
         cout << "Only coded for Lagrange polynomial approximation\n";
         cout << key << endl;
         throw;
      }
   }

   // Face
   if (enti->getLevel() == 2)
   {
      xApproxFunctionScalarTriLagrange* shapeFuncF = dynamic_cast<xApproxFunctionScalarTriLagrange*>((func.get()));
      if (shapeFuncF)
      {
         std::vector<double> nodeuvw;
         shapeFuncF->getNodeUvw(nodeuvw);
         return xtensor::xPoint(nodeuvw[0], nodeuvw[1], 0);
      }
      else
      {
         cout << "Only coded for Lagrange polynomial approximation\n";
         cout << key << endl;
         throw;
      }
   }
   else
   {
      cout << "xSpaceXFEMShiftedCorrected::getInterpolantUVW Only Coded for entity level 0, 1, 2 " << __FILE__ << __LINE__
           << std::endl;
      throw;
   }
}

}  // namespace xfem
