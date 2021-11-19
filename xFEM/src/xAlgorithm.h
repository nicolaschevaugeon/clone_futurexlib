/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _FORMULATION_H
#define _FORMULATION_H

#include <algorithm>
#include <ctime>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "xApplyCommandOnIntegrationRule.h"
#include "xAssembler.h"
#include "xCommandOnGeomElem.h"
#include "xDebug.h"
#include "xEntityToEntity.h"
#include "xEnv.h"
#include "xEval.h"
#include "xField.h"
#include "xFiniteElement.h"
#include "xForm.h"
#include "xGeomElem.h"
#include "xIntegrationRule.h"
#include "xOperators.h"
#include "xRegion.h"
#include "xSolVisitor.h"
#include "xStateOfValue.h"
#include "xValue.h"
#include "xValueManager.h"

// AOMD
#include "xValueLinearCombination.h"

namespace xfem
{
// defined algorithm :
// DeclareInterpolation
// DeclareInterpolationAndStateMultLag
// DeleteInterpolation
// UpdateInterpolation
// DeclareValueField
// DeleteValueField
// ApplyCommandOnIntegrationRule
// Visit
// Accumulate
// Assemble
// DirichletBoundaryCondition
// DeclareCylindricalDirichletBoundaryCondition3D
// DeclareState
// DeleteState
// L2Projection
// generateMeshToMeshConstrainWithL2Projection
// DeclareInterpolationHanging
// FillField
// ComputeNnz
// AssembleGraph
// ComputResultant
// ComputMoment

/*!
 *
 * Needed
 * FIELD member :
 *     xValueManagerDist<double> * getValueManager()
 *     std::vector<spacePtr>::iterator begin()
 *     std::vector<spacePtr>::iterator end()
 *
 * CREA :
 *     could be a object function, a boost or std function, a lambda function, ....
 *     In all cases must respect following signature :
 *         xValue<typename FIELD::value_t>* (const xValKey&) const
 */
template <class FIELD, class CREA, class ITER>
void DeclareInterpolation(const FIELD& fct, const CREA& value_creator, ITER it, ITER end)
{
   const bool debug = xdebug_flag;
   if (debug) std::cout << "begin declare\n";
   std::set<xValKey> remain;
   auto double_manager = fct.getValueManager();
   //    ITER it_save = it;
   for (; it != end; ++it)
   {
      if (debug) std::cout << "elt in  declare\n";
      AOMD::mEntity* e = *it;
      xFiniteElement FEM;
      FEM.setKeys(e, fct.begin(), fct.end());
      if (debug) std::cout << "size of fct " << fct.size() << std::endl;
      if (debug) std::cout << "size of key " << FEM.sizeKey() << std::endl;
      for (xFiniteElement::iterKey itk = FEM.beginKey(); itk != FEM.endKey(); ++itk)
      {
         if (!double_manager->insert(*itk, value_creator)) remain.insert(*itk);
      }
   }
   if (debug) std::cout << "remain size " << remain.size() << std::endl;
   while (!remain.empty())
   {
      // see effective stl of Scott Meyers page 45 to understand the line below
      for (std::set<xValKey>::iterator itk = remain.begin(); itk != remain.end();)
      {
         if (debug) std::cout << *itk << std::endl;
         if (double_manager->insert(*itk, value_creator))
         {
            if (debug) std::cout << "erasing " << std::endl;
            remain.erase(itk++);
            if (debug) std::cout << "after erasing " << std::endl;
         }
         else
            ++itk;
      }
      std::cout << "size of remain " << remain.size() << std::endl;
   }
   if (debug) std::cout << "end declare\n";

   //  ITER ite = it_save;
   // for(; ite != end; ++ite)
   //  {
   //    AOMD::mEntity*e  = *ite ;
   //    fct.storeApproximation(e);
   //  }

   return;
}

template <class FIELD, class CREA, class ITER>
void DeclareInterpolation2(const FIELD& fct, const CREA& value_creator, ITER it, ITER end)
{
   std::set<xValKey> remain;
   auto double_manager = fct.getValueManager();
   std::vector<xValKey> glob;
   std::vector<xValKey> loc;
   auto space_range = xtool::make_range(fct.begin(), fct.end());
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      glob.clear();
      for (auto space : space_range)
      {
         loc.clear();
         space->getKeys(e, &loc);
         glob.insert(glob.end(), loc.begin(), loc.end());
      }
      for (auto& key : glob)
      {
         if (!double_manager->insert(key, value_creator)) remain.insert(key);
      }
   }
   while (!remain.empty())
   {
      // see effective stl of Scott Meyers page 45 to understand the line below
      for (std::set<xValKey>::iterator itk = remain.begin(); itk != remain.end();)
      {
         if (double_manager->insert(*itk, value_creator))
         {
            remain.erase(itk++);
         }
         else
            ++itk;
      }
   }

   return;
}
// function deprecated
// use the classical "DeclareInterpolation" and "DeclareState" framework with a xValueCreatorLinkOnFront creator
/*template <class FIELD>
void DeclareInterpolationAndStateMultLag (const FIELD& fct,
        const xMesh* bnd,
        std::function<xValue<double>* (const xValKey&, xValue<double>*)> set_dof)
{
    const bool debug = xdebug_flag;
    xEntityToEntity interf2appro = xCreator();

    xValueManagerDist<double>* double_manager = fct.getValueManager();

    if (bnd->size(0) == 0) return;

    AOMD::mEntity* e_dummy = *bnd->begin(0);
    xFiniteElement FEM;
    FEM.setKeys(e_dummy, fct.begin(), fct.end());
    xValueCreator<xValueDouble> value_creator;

    //set of nodes already visited
    mMeshEntityContainer visited_nodes;

    //first loop to select the vital edge
    for(xIter it = bnd->begin(0); it != bnd->end(0); ++it)
    {
        AOMD::mEntity * e_bnd = *it;
        AOMD::mEntity * e = interf2appro(e_bnd);

        if (e->getLevel() == 0)
        {
            // create a dof
            std::for_each(FEM.beginKey(), FEM.endKey(),bind2nd(mem_fun_ref(&xValKey::setEnti), e));
            std::vector<xValue<double>*> vals;
            double_manager->insert(FEM.beginKey(), FEM.endKey(), value_creator, vals);
            std::transform(FEM.beginKey(), FEM.endKey(), vals.begin(), vals.begin(), set_dof);
            visited_nodes.add(e);
        }
        else if(e->getLevel() == 1)
        {
            //check if 0, 1 or 2 nodes have dofs
            //if 0, create two dofs and link them
            AOMD::mEntity* v1 = e->get(0,0);
            AOMD::mEntity* v2 = e->get(0,1);
            if (!visited_nodes.find(v1) && !visited_nodes.find(v2))
            {
                std::for_each(FEM.beginKey(), FEM.endKey(),bind2nd(mem_fun_ref(&xValKey::setEnti), v1));
                std::vector<xValue<double>*> vals;
                double_manager->insert(FEM.beginKey(), FEM.endKey(), value_creator, vals);
                std::transform(FEM.beginKey(), FEM.endKey(), vals.begin(), vals.begin(), set_dof);

                std::for_each(FEM.beginKey(), FEM.endKey(),bind2nd(mem_fun_ref(&xValKey::setEnti), v2));
                std::vector<xValue<double>*>::const_iterator itv =  vals.begin();
                for (xFiniteElement::iterKey itk = FEM.beginKey(); itk != FEM.endKey(); ++itk, ++itv)
                {
                    xValue<double>* val2 = new xValueLinearCombination(1., *itv);
                    double_manager->insert(*itk, val2);
                }
                visited_nodes.add(v1);
                visited_nodes.add(v2);
            }
        }
    }
    //second loop to set up the final dofs.
    for(xIter it = bnd->begin(0); it != bnd->end(0); ++it)
    {
        AOMD::mEntity * e_bnd = *it;
        AOMD::mEntity * e = interf2appro(e_bnd);
        if(e->getLevel() == 1)
        {
            //check if 0, 1 or 2 nodes have dofs
            //if 1, create the other dof and link it
            AOMD::mEntity* v1 = e->get(0,0);
            AOMD::mEntity* v2 = e->get(0,1);
            bool b1 = (visited_nodes.find(v1));
            bool b2 = (visited_nodes.find(v2));
            assert(b1 || b2);
            if (!b1 || !b2)
            {
                if (b2) swap(v1, v2);
                std::for_each(FEM.beginKey(), FEM.endKey(),bind2nd(mem_fun_ref(&xValKey::setEnti), v1));
                std::vector<xValue<double>*> vals;
                double_manager->getValPtr(FEM.beginKey(), FEM.endKey(), vals);

                std::for_each(FEM.beginKey(), FEM.endKey(),bind2nd(mem_fun_ref(&xValKey::setEnti), v2));
                std::vector<xValue<double>*>::const_iterator itv = vals.begin();
                for (xFiniteElement::iterKey itk = FEM.beginKey(); itk != FEM.endKey(); ++itk, ++itv)
                {
                    xValue<double>* val2 = new xValueLinearCombination(1., *itv);
                    double_manager->insert(*itk, val2);
                }
                visited_nodes.add(v2);
            }
        }
    }
    return;
 };*/

template <class FIELD, class ITER>
void DeleteInterpolation(const FIELD& fct, ITER it, ITER end)
{
   auto double_manager = fct.getValueManager();
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      xFiniteElement FEM;
      FEM.setKeys(e, fct.begin(), fct.end());
      for (xFiniteElement::iterKey itk = FEM.beginKey(); itk != FEM.endKey(); ++itk) double_manager->erase(*itk);
   }
   // fct.deleteApproximation();
}

/// This function update a already declared interpolation with the use of a dedicated updator.
//! Be careful to call this function with a std::ref(value_updator) instead of just value_updator. This
//! prevent the copy of your object value_updator and the use of this copy in this function. In many case it
//! doesn't matter to use the object or its copy but if you wish to use the object after this function you must
//! use this reference mecanism.
template <class FIELD, class UPDATE, class ITER>
void UpdateInterpolation(const FIELD& fct, const UPDATE& value_updator, ITER it, ITER end)
{
   const bool debug = xdebug_flag;
   if (debug) std::cout << "begin interpolation update\n";
   std::set<xValKey> remain;
   auto double_manager = fct.getValueManager();
   //    ITER it_save = it;
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      xFiniteElement FEM;
      FEM.setKeys(e, fct.begin(), fct.end());
      for (xFiniteElement::iterKey itk = FEM.beginKey(); itk != FEM.endKey(); ++itk)
      {
         if (!double_manager->update(*itk, value_updator)) remain.insert(*itk);
      }
   }
   while (!remain.empty())
   {
      // see effective stl of Scott Meyers page 45 to understand the line below
      for (std::set<xValKey>::iterator itk = remain.begin(); itk != remain.end();)
      {
         if (debug) std::cout << *itk << std::endl;
         if (double_manager->update(*itk, value_updator))
         {
            if (debug) std::cout << "erasing " << std::endl;
            remain.erase(itk++);
            if (debug) std::cout << "after erasing " << std::endl;
         }
         else
            ++itk;
      }
      std::cout << "size of remain " << remain.size() << std::endl;
   }
   if (debug) std::cout << "end declare\n";

   return;
}

// ?? does this function need to exist ???
// now that DeclareInterpolation is template on value type ....
template <class Field, class CREA, class ITER>
void DeclareValueField(const Field& fct, const CREA& value_creator, ITER it, ITER end)
{
   const bool debug = xdebug_flag;
   std::set<xValKey> remain;
   typename Field::value_manager_t* val_manager = fct.getValueManager();

   // ITER it_save = it;

   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      xFiniteElementKeysOnly FEM;
      FEM.setKeys(e, fct.begin(), fct.end());
      for (xFiniteElementKeysOnly::iterKey itk = FEM.beginKey(); itk != FEM.endKey(); ++itk)
      {
         if (!val_manager->insert(*itk, value_creator)) remain.insert(*itk);
      }
   }
   while (!remain.empty())
   {
      // see effective stl of Scott Meyers page 45 to understand the line below
      for (std::set<xValKey>::iterator itk = remain.begin(); itk != remain.end();)
      {
         if (val_manager->insert(*itk, value_creator))
         {
            if (debug) std::cout << "erasing " << std::endl;
            remain.erase(itk++);
            if (debug) std::cout << "after erasing " << std::endl;
         }
         else
            ++itk;
      }
      if (debug) std::cout << "size of remain " << remain.size() << std::endl;
   }
}

template <class Field, class ITER>
void DeleteValueField(const Field& fct, ITER it, ITER end)
{
   typename Field::value_manager_t* val_manager = fct.getValueManager();
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      xFiniteElementKeysOnly FEM;
      FEM.setKeys(e, fct.begin(), fct.end());
      for (xFiniteElementKeysOnly::iterKey itk = FEM.beginKey(); itk != FEM.endKey(); ++itk) val_manager->erase(*itk);
   }
}

template <class V, class ITER>
void Visit(V& visitor, ITER first, ITER last)
{
   std::transform(first, last, first, visitor);
}

template <class V, class ITER>
void Visit(const V& visitor, ITER first, ITER last)
{
   std::transform(first, last, first, visitor);
}

template <class BinaryOperator, class ITER, class T>
void Accumulate(BinaryOperator& accumulator, ITER first, ITER last, T& res)
{
   res = std::accumulate(first, last, res, accumulator);
}

template <class BinaryOperator, class ITER, class T>
void Accumulate(const BinaryOperator& accumulator, ITER first, ITER last, T& res)
{
   res = std::accumulate(first, last, res, accumulator);
}

template <class BILINEARFORM, class ASSEMBLER, class FIELD, class ITER>
void Assemble(BILINEARFORM& bilin, ASSEMBLER& assembler, const xIntegrationRule& integration_rule, const FIELD& test,
              const FIELD& trial, ITER it, ITER end, xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   Assemble(bilin, assembler, integration_rule, test, trial, it, end, integ2appro, integ2appro);
}

template <class BILINEARFORM, class ASSEMBLER, class FIELD, class ITER>
void Assemble(BILINEARFORM& bilin, ASSEMBLER& assembler, const xIntegrationRule& integration_rule, const FIELD& test,
              const FIELD& trial, ITER it, ITER end, xEntityToEntity integ2appro_test, xEntityToEntity integ2appro_trial)
{
   typedef femFcts_t femFcts;
   const bool debug = false;  // xDebugSingleton::instance().flag;
   xIntegrateFormCommand command(&bilin);
   if (&trial == &test)
   {
      if (debug) std::cout << "Assemble : TEST==TRIAL\n";
      for (; it != end; ++it)
      {
         AOMD::mEntity* e_integ = *it;
         AOMD::mEntity* e_appro = integ2appro_test(e_integ);
         femFcts trialfcts(trial.beginFcts(e_appro), trial.endFcts(e_appro));
         if (!trialfcts.empty())
         {
            if (debug) std::cout << "assembler.ass-1\n";
            xGeomElem geo_appro(e_appro);
            bilin.init(&geo_appro, &trialfcts);
            integration_rule.accept(command, e_integ);
            bilin.finalize();
            if (debug) std::cout << "assembler.assemble\n";
            assembler.assemble_matrix(test.beginValues(e_appro), test.endValues(e_appro), bilin.getMatrix());
         }
      }
   }
   else
   {
      if (debug) std::cout << "Assemble : TEST # TRIAL\n";
      for (; it != end; ++it)
      {
         AOMD::mEntity* e_integ = *it;
         AOMD::mEntity* e_appro_test = integ2appro_test(e_integ);
         femFcts testfcts(test.beginFcts(e_appro_test), test.endFcts(e_appro_test));
         AOMD::mEntity* e_appro_trial = integ2appro_trial(e_integ);
         femFcts trialfcts(trial.beginFcts(e_appro_trial), trial.endFcts(e_appro_trial));
         if ((!testfcts.empty()) && !(trialfcts.empty()))
         {
            xGeomElem geo_appro_test(e_appro_test);
            xGeomElem geo_appro_trial(e_appro_trial);
            bilin.init(&geo_appro_test, &testfcts, &geo_appro_trial, &trialfcts);
            integration_rule.accept(command, e_integ);
            bilin.finalize();
            assembler.assemble_matrix(test.beginValues(e_appro_test), test.endValues(e_appro_test),
                                      trial.beginValues(e_appro_trial), trial.endValues(e_appro_trial), bilin.getMatrix());
         }
      }
   }
   return;
}  // end Assemble xFormBilinear

template <class LINEARFORM, class ASSEMBLER, class FIELD, class ITER>
void Assemble(LINEARFORM& lin, ASSEMBLER& assembler, const xIntegrationRule& integration_rule, const FIELD& test, ITER it,
              ITER end, xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   const bool debug = false;
   xIntegrateFormCommand command(&lin);
   for (; it != end; ++it)
   {
      AOMD::mEntity* e_integ = *it;
      if (debug)
      {
         std::cout << "In linear assemble" << std::endl;
         std::cout << "e_integ is " << std::endl;
         e_integ->print();
      }
      AOMD::mEntity* e_appro = integ2appro(e_integ);
      if (debug)
      {
         std::cout << "e_appro is " << std::endl;
         e_appro->print();
      }
      femFcts_t fcts(test.beginFcts(e_appro), test.endFcts(e_appro));
      if (!fcts.empty())
      {
         xGeomElem geo_appro(e_appro);
         command.openApproxElem(&geo_appro);
         lin.init(&geo_appro, &fcts);
         integration_rule.accept(command, e_integ);
         assembler.assemble_vector(test.beginValues(e_appro), test.endValues(e_appro), lin.getVector());
         lin.finalize();
         command.closeApproxElem(&geo_appro);
         if (debug) std::cout << "assembling " << lin.getVector() << std::endl;
      }
   }
   return;
}

template <class ZEROFORM, class ASSEMBLER, class ITER>
void Assemble(ZEROFORM& zer, ASSEMBLER& assembler, const xIntegrationRule& integration_rule, ITER it, ITER end,
              xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   xIntegrateFormCommand command(&zer);
   for (; it != end; ++it)
   {
      AOMD::mEntity* e_integ = *it;
      AOMD::mEntity* e_appro = integ2appro(e_integ);
      xGeomElem geo_appro(e_appro);
      zer.init(&geo_appro);
      integration_rule.accept(command, e_integ);
      zer.finalize();
      assembler.assemble_scalar(zer.getScalar().getVal());
   }
   return;
}

template <class ZEROFORM, class ITER>
void Assemble(ZEROFORM& zer, double& result, const xIntegrationRule& integration_rule, ITER it, ITER end,
              xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   xIntegrateFormCommand command(&zer);
   for (; it != end; ++it)
   {
      AOMD::mEntity* e_integ = *it;
      AOMD::mEntity* e_appro = integ2appro(e_integ);
      xGeomElem geo_appro(e_appro);
      zer.init(&geo_appro);
      integration_rule.accept(command, e_integ);
      zer.finalize();
      result += zer.getScalar().getVal();  // hans
   }
   return;
}

template <class FIELD, class ITER>
void DirichletBoundaryCondition(const FIELD& listFunctionSpaceconn, const string& phys, ITER it, ITER end,
                                typename FIELD::value_t val = xtool::xDataType<typename FIELD::value_t>::zero())
{
   const bool debug = xdebug_flag;
   if (debug) std::cout << "in  DirichletBoundaryCondition " << std::endl;
   if (debug) std::cout << "phys " << phys << std::endl;
   auto double_manager = listFunctionSpaceconn.getValueManager();
   xIsPhys f_phys(phys);

   for (; it != end; ++it)
   {
      if (debug) std::cout << "element in DirichletBoundaryCondition" << std::endl;
      AOMD::mEntity* e = *it;
      xFiniteElement FEM;
      FEM.setKeys(e, listFunctionSpaceconn.begin(), listFunctionSpaceconn.end());
      if (debug) std::cout << "size of keys before filter " << FEM.sizeKey() << std::endl;
      FEM.filterKey(f_phys);

      typename FIELD::value_container_t vals;
      double_manager->getValPtr(FEM.beginKey(), FEM.endKey(), vals);
      if (debug) std::cout << "size of keys after  filter " << FEM.sizeKey() << std::endl;
      // set the value
      for (auto it = vals.begin(); it != vals.end(); ++it)
      {
         (*it)->setVal(val);
      }
      std::transform(FEM.beginKey(), FEM.endKey(), vals.begin(), vals.begin(), xStateFixedCreator<typename FIELD::value_t>());

      /* The following is an unsuccessful attempt to set dirichlet bc when more than one dof define the value of
           the current field at one entity...
           if (vals.size()==1) {
           vals[0]->setVal(val);
           vals[0]->setState(new xStateOfValueFixed(vals[0]));
      //      std::transform(FEM.beginKey(), FEM.endKey(), vals.begin(), vals.begin(),
      //		   xStateFixedCreator());
      }
      else if (vals.size()>1) {
      femFcts_t functions =  FEM.getFcts();
      std::vector<double > coefs(vals.size()-1);
      std::fill(coefs.begin(),coefs.end(),-1.);
      std::vector<xValue<double>*> vals2(vals.size()-1);
      std::copy((++vals.begin()), vals.end(),vals2.begin());
      xValueLinearCombination * comb=  new xValueLinearCombination( coefs,  vals2, val);
      vals[0]->setState(new xStateOfValueLinearCombination(comb) );
      }*/
   }
}

/*!
 *
 * Needed
 * FIELD member :
 *     xValueManagerDist<double> * getValueManager()
 *     std::vector<spacePtr>::iterator begin()
 *     std::vector<spacePtr>::iterator end()
 * CREA :
 *     could be a object function, a boost or std function, a lambda function, ....
 *     In all cases must respect following signature :
 *         xValue<typename FIELD::value_t>* (const xValKey&) const
 */
template <class FIELD, class CREA, class ITER>
void DeclareInterpolationRigid(const FIELD& listFunctionSpaceconn, const CREA& value_creator, const string& phys, ITER it,
                               ITER end)
{
   const bool debug = xdebug_flag;
   if (debug) std::cout << "in  DirichletBoundaryCondition " << std::endl;
   if (debug) std::cout << "phys " << phys << std::endl;
   auto double_manager = listFunctionSpaceconn.getValueManager();

   bool first = true;
   xValue<typename FIELD::value_t>* val_dof = nullptr;
   for (; it != end; ++it)
   {
      if (debug) std::cout << "element in DirichletBoundaryCondition" << std::endl;
      AOMD::mEntity* e = *it;
      xFiniteElement FEM;
      FEM.setKeys(e, listFunctionSpaceconn.begin(), listFunctionSpaceconn.end());
      if (debug) std::cout << "size of keys before filter " << FEM.sizeKey() << std::endl;
      FEM.filterKey(xIsPhys(phys));

      if (debug) std::cout << "size of keys after  filter " << FEM.sizeKey() << std::endl;
      // set the value
      xFiniteElement::iterKey itk = FEM.beginKey();
      xFiniteElement::iterKey ite = FEM.endKey();

      for (; itk != ite; ++itk)
      {
         const xValKey& key = *itk;
         if (first)
         {
            val_dof = double_manager->insert(key, value_creator);
            assert(val_dof);
            first = false;
         }
         else
         {
            assert(val_dof);
            double_manager->insert(key, new xValueLinearCombination<typename FIELD::value_t>(
                                            xtool::xDataType<typename FIELD::value_t>::one(), val_dof));
         }
      }
   }
}

template <class FIELD, class ITER>
void TreatmentOfEssentialEnvCylindrical(const FIELD& fct,
                                        std::function<xValue<typename FIELD::value_t>*(const xValKey&)> value_creator, ITER itenv,
                                        ITER endenv, const std::vector<string>& disp_3D, const std::vector<string>& vars,
                                        const xtensor::xPoint& origin_coord, const xtensor::xVector<>& z_direction, xMesh* mesh)
{
   // const bool debug = xdebug_flag;
   xValKey::ids_size_t phys_r = xKeyInfo::getPhysId("DISPLACEMENT_R");
   xValKey::ids_size_t phys_t = xKeyInfo::getPhysId("DISPLACEMENT_TETA");
   for (; itenv != endenv; ++itenv)
   {
      const xEnv& env = *itenv;
      string phys = env.Phys;
      int type = env.Type;
      int entity = env.Entity;
      if (find(vars.begin(), vars.end(), phys) != vars.end())
      {
         if (type == FIX)
         {
            xValKey::ids_size_t physid = xKeyInfo::getPhysId(phys);
            xClassRegion bc(mesh, entity, env.getDimension());
            double val = env.getValue();
            DeclareCylindricalDirichletBoundaryCondition3D(fct, value_creator, origin_coord, z_direction,
                                                           (phys_r == physid) ? 0 : ((phys_t == physid) ? 1 : 2), disp_3D[0],
                                                           disp_3D[1], disp_3D[2], bc.begin(), bc.end(), val);
         }
      }
   }
   return;
}

/*!
 * Giving a point and a z direction describing partially a Cylindrical coordinate
 * this function chose an arbitrary Cylindrical coordinate (by choosing r direction)
 * and impose in this coordinate Dirichlet type boundary condition.
 * Needed
 * FIELD member :
 *     xValueManagerDist<double> * getValueManager()
 *     std::vector<spacePtr>::iterator begin()
 *     std::vector<spacePtr>::iterator end()
 *     type : FIELD::value_t
 * CREA :
 *     could be a object function, a boost or std function, a lambda function, ....
 *     In all cases must respect following signature :
 *         xValue<typename FIELD::value_t>* (const xValKey&) const
 */
template <class FIELD, class CREA, class ITER>
void DeclareCylindricalDirichletBoundaryCondition3D(const FIELD& field, const CREA& value_creator,
                                                    const xtensor::xPoint& origin_coord, xtensor::xVector<> z_direction,
                                                    const int fixed_componant, const string& phys_x, const string& phys_y,
                                                    const string& phys_z, ITER it, ITER end,
                                                    double val = xtool::xDataType<typename FIELD::value_t>::zero())
{
   auto double_manager = field.getValueManager();

   // instantiate predicate for phys filtering
   xIsPhys f_phys_x(phys_x);
   xIsPhys f_phys_y(phys_y);
   xIsPhys f_phys_z(phys_z);

   // normalize z_direction
   z_direction.norm();

   // component of z_direction
   const double v0 = z_direction(0);
   const double v1 = z_direction(1);
   const double v2 = z_direction(2);

   // z_direction component at power 2
   const double v02 = v0 * v0;
   const double v12 = v1 * v1;
   const double v22 = v2 * v2;

   // projection of z_direction on global plane yz
   double nt1 = v12 + v22;
   // projection of z_direction on global plane xz
   double nt2 = v02 + v22;
   // projection of z_direction on global plane xy
   double nt3 = v02 + v12;

   // constructing pseudo_x_direction using biggest z_direction projection on planes to construct pseudo orthogonal vector
   xtensor::xVector<> pseudo_x_direction;
   if (nt1 > nt3)
   {
      if (nt1 > nt2)
      {
         nt1 = sqrt(nt1);
         pseudo_x_direction(1) = v2 / nt1;
         pseudo_x_direction(2) = -v1 / nt1;
      }
      else
      {
         nt2 = sqrt(nt2);
         pseudo_x_direction(0) = -v2 / nt2;
         pseudo_x_direction(2) = v0 / nt2;
      }
   }
   else
   {
      if (nt3 > nt2)
      {
         nt3 = sqrt(nt3);
         pseudo_x_direction(0) = v1 / nt3;
         pseudo_x_direction(1) = -v0 / nt3;
      }
      else
      {
         nt2 = sqrt(nt2);
         pseudo_x_direction(0) = -v2 / nt2;
         pseudo_x_direction(2) = v0 / nt2;
      }
   }

   // construct y_direction as orthogonal vector to z_direction and pseudo_x_direction
   xtensor::xVector<> y_direction = z_direction % pseudo_x_direction;
   y_direction.norm();

   // finalizing construction by computing x_direction as orthogonal vector to y_direction and z_direction
   xtensor::xVector<> x_direction = y_direction % z_direction;
   x_direction.norm();

   // x_direction,y_direction and z_direction correspond to base of user Cartesian coordinate (having z and origin corresponding
   // to cylindrical coordinate). They are used also as component of the rotation matrix to transform displacement UGN of a point
   // N in global coordinate into UN , displacement of point N in user coordinate. UNx=x_direction.UGN UNy=y_direction.UGN
   // UNz=z_direction.UGN
   //

   xStateFixedCreator<typename FIELD::value_t> fixed_value_creator;
   int i, nbkey, col_id, slvb, slv1, slv2;
   double epsilon = std::numeric_limits<typename FIELD::value_t>::epsilon();
   double x_cart_cyl, y_cart_cyl;
   // double z_cart_cyl;
   double alpha, beta;
   // double gama0,gama1,gama2
   double cst;
   double *x_max, *y_max;
   double pivot, csteq;
   typename xValueLinearCombination<typename FIELD::value_t>::coeffs_t coeffs;
   typename xValueLinearCombination<typename FIELD::value_t>::xvalues_t master;
   xtensor::xVector<> vgama;
   xFiniteElementKeysOnly FEM;
   xFiniteElementKeysOnly::iterKey itk;
   xFiniteElementKeysOnly::iterKey itke;
   xValKey keyx, keyy, keyz, keys, key_coef[2];

   // loop on element where cylindrical Dirichlet have to be set
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      // Only x,y,z triplet considered for nodes as we only have coordinates of this kind of entity
      // Loops on them
      const int nbn = e->size(0);
      for (i = 0; i < nbn; ++i)
      {
         AOMD::mEntity* ne = e->get(0, i);

         // get keys
         FEM.setKeys(ne, field.begin(), field.end());
         nbkey = FEM.sizeKey();
         // at least one key for this node, check if value created or not
         if (nbkey)
         {
            // if at least one value exist, skip creation for this node
            bool skip = false;
            for (itk = FEM.beginKey(), itke = FEM.endKey(); itk != itke; ++itk)
               if (double_manager->find(*itk))
               {
                  skip = true;
                  break;
               }

            if (skip) continue;
         }
         // no key for this node .... skip it
         else
            continue;

         // retrieve global x,y,z corresponding key (this is done as field may for example be enriched
         // and we only treat here normal dof : TODO add enrichment )
         for (itk = FEM.beginKey(), itke = FEM.endKey(); itk != itke; ++itk)
         {
            xValKey& key = (*itk);
            if (f_phys_x(key))
               keyx = key;
            else if (f_phys_y(key))
               keyy = key;
            else if (f_phys_z(key))
               keyz = key;
         }

         // compute N coordinate in Cartesian user coordinate
         xtensor::xPoint N_point = (static_cast<AOMD::mVertex*>(ne))->point();
         xtensor::xVector<> N_vect_from_cyl_orig(origin_coord, N_point);
         x_cart_cyl = x_direction * N_vect_from_cyl_orig;
         y_cart_cyl = y_direction * N_vect_from_cyl_orig;
         // z_cart_cyl=z_direction*N_vect_from_cyl_orig;
         bool y_not_null = (y_cart_cyl > epsilon || y_cart_cyl < -epsilon);
         bool x_not_null = (x_cart_cyl > epsilon || x_cart_cyl < -epsilon);

         // reset
         coeffs.resize(2);
         master.resize(2);
         x_max = y_max = nullptr;
         slvb = slv1 = slv2 = -1;
         pivot = 0.;
         csteq = 0.;
         col_id = 0;

         switch (fixed_componant)
         {
            // r component blocked
            case 0:
            {
               if (y_not_null || x_not_null)
               {
                  if (y_not_null)
                  {
                     // compute alpha,beta of the linear relation UNy = alpha + beta * UNx where UNx and UNy are
                     // displacement of N in user Cartesian coordinate.
                     beta = 1. / y_cart_cyl;
                     if (val > 0.)
                     {
                        alpha = val * sqrt(x_cart_cyl * x_cart_cyl + y_cart_cyl * y_cart_cyl) * beta;
                     }
                     else
                        alpha = 0.0;
                     beta *= -x_cart_cyl;

                     // compute linear relation in global Cartesian coordinate by applying rotation from global
                     // to user Cartesian coordinate to UGN and replace in linear equation above.
                     // Rewriting equation above this way UNy-beta*UNx=alpha
                     // we have
                     //  y_direction.UGN-x_direction.UGN*beta=alpha
                     //  (y_direction-x_direction.beta).UGN=alpha
                     vgama = y_direction - (x_direction * beta);
                  }
                  else
                  {
                     // y_cart_cyl is null but as x_cart_cyl is not null UNy and
                     // UNz are free but UNx=alpha => x_direction.UGN=alpha
                     vgama = x_direction;
                     if (x_cart_cyl < 0.)
                        alpha = -val;
                     else
                        alpha = val;
                  }
               }
               // if x_cart_cyl and y_cart_cyl null :
               //   UNz is free but UNx=alpha=val and UNy=0
               //
               //   =>
               //                x_direction.UGN=alpha
               //                y_direction.UGN=0
               //
               //  This form a system of 2 equations with 3 unknowns which impose 2 linear combination
               //  making 2 unknowns slave of the 3th one.
               //  It can be fully computed (i.e. both 2 unknows depend only on the 3th one) or one unknows depend only of the 3th
               //  one and the other unknow depens on other 2. This is this last solution which have been chossen leading to
               //  relation depending on other
               else
               {
                  // search first pivot (max of all direction vector coefficiant)
                  x_max = std::max_element(&x_direction(0), (&x_direction(0)) + 3);
                  y_max = std::max_element(&y_direction(0), (&y_direction(0)) + 3);

                  if (*x_max > *y_max)
                  {
                     col_id = x_max - (&x_direction(0));
                     assert(col_id > -1 && col_id < 3);
                     csteq = 1. / (*x_max);
                     pivot = -y_direction(col_id) * csteq;
                     vgama = y_direction + x_direction * pivot;
                     alpha = val * pivot;
                  }
                  else
                  {
                     col_id = y_max - (&y_direction(0));
                     assert(col_id > -1 && col_id < 3);
                     csteq = 1. / (*y_max);
                     pivot = -x_direction(col_id) * csteq;
                     vgama = x_direction + y_direction * pivot;
                     alpha = val;
                  }
                  int id = 0;
                  // fixing slave dof of first equation used as pivot
                  while (slv2 < 0)
                  {
                     assert(id < 3);
                     if (id != col_id)
                     {
                        if (slv1 < 0)
                        {
                           slv1 = id;
                        }
                        else
                           slv2 = id;
                        id++;
                     }
                  }
               }
               break;
            }
            // teta composant blocked
            case 1:
            {
               if (y_not_null)
               {
                  // compute alpha,beta of the linear relation UNx = alpha + beta * UNy where UNx and UNy are
                  // displacement of N in user Cartesian coordinate.
                  beta = 1. / y_cart_cyl;
                  if (val > 0.)
                  {
                     alpha = -val * (x_cart_cyl * x_cart_cyl + y_cart_cyl * y_cart_cyl) * beta;
                  }
                  else
                     alpha = 0.0;
                  beta *= x_cart_cyl;

                  // compute linear relation in global Cartesian coordinate by applying rotation from global
                  // to user Cartesian coordinate to UGN and replace in linear equation above.
                  // Rewriting equation above this way UNx-beta*UNy=alpha
                  // we have
                  //  x_direction.UGN-y_direction.UGN*beta=alpha
                  //  (x_direction-y_direction.beta).UGN=alpha
                  vgama = x_direction - (y_direction * beta);
               }
               else
               {
                  // y_cart_cyl is null
                  // UNx and UNz are free but UNy=alpha => y_direction.UGN=alpha
                  // if x_cart_cyl is not null  : alpha=+/- r*teta imposed
                  // if x_cart_cyl is null : apha=0.
                  vgama = y_direction;
                  if (x_not_null)
                  {
                     alpha = val * sqrt(x_cart_cyl * x_cart_cyl + y_cart_cyl * y_cart_cyl);
                     if (x_cart_cyl < 0.) alpha = -alpha;
                  }
                  else
                     alpha = 0.;
               }
               break;
            }
            // z composant blocked
            case 2:
            {
               // UNx and UNz are free but UNz=alpha => z_direction.UGN=alpha=val
               vgama = z_direction;
               alpha = val;
               break;
            }
            // default to avoid -Wmaybe-uninitialized
            default:
               throw -2354;
         }

         // One component have to be eliminated : choosing biggest vgama coefficient for slave dof
         const double gama0 = fabs(vgama(0));
         const double gama1 = fabs(vgama(1));
         const double gama2 = fabs(vgama(2));
         if (gama0 > gama1)
         {
            // x biggest
            if (gama0 > gama2)
            {
               slvb = 0;
               keys = keyx;
               key_coef[0] = keyy;
               key_coef[1] = keyz;
               cst = 1. / vgama(0);
               coeffs[0] = -vgama(1) * cst;
               coeffs[1] = -vgama(2) * cst;
            }
            // z biggest
            else
            {
               slvb = 2;
               keys = keyz;
               key_coef[0] = keyx;
               key_coef[1] = keyy;
               cst = 1. / vgama(2);
               coeffs[0] = -vgama(0) * cst;
               coeffs[1] = -vgama(1) * cst;
            }
         }
         else
         {
            // y biggest
            if (gama1 > gama2)
            {
               slvb = 1;
               keys = keyy;
               key_coef[0] = keyx;
               key_coef[1] = keyz;
               cst = 1. / vgama(1);
               coeffs[0] = -vgama(0) * cst;
               coeffs[1] = -vgama(2) * cst;
            }
            // z biggest
            else
            {
               slvb = 2;
               keys = keyz;
               key_coef[0] = keyx;
               key_coef[1] = keyy;
               cst = 1. / vgama(2);
               coeffs[0] = -vgama(0) * cst;
               coeffs[1] = -vgama(1) * cst;
            }
         }

         // scale cst
         cst *= alpha;

         // filter null coefficient and creat free value
         // if second coef is null
         if (fabs(coeffs[1]) < epsilon)
         {
            // if both are null just resize to zero
            if (fabs(coeffs[0]) < epsilon)
            {
               coeffs.resize(0);
            }
            // otherwise remove last by resize to 1 and create value
            else
            {
               coeffs.resize(1);
               master.resize(1);
               master[0] = double_manager->insert(key_coef[0], value_creator);
            }
         }
         else
         {
            // if first is null shift, resize and create value
            if (fabs(coeffs[0]) < epsilon)
            {
               coeffs[0] = coeffs[1];
               coeffs.resize(1);
               master.resize(1);
               master[0] = double_manager->insert(key_coef[1], value_creator);
            }
            // if no null coef create 2 values
            else
            {
               master[0] = double_manager->insert(key_coef[0], value_creator);
               master[1] = double_manager->insert(key_coef[1], value_creator);
            }
         }
         xValue<double>* vs;

         // depending if linear relation reduce to constant or not generate last value
         if (coeffs.size())
         {
            // create slave value
            vs = static_cast<xValue<double>*>(new xValueLinearCombination<double>(coeffs, master, cst));
            double_manager->insert(keys, vs);
         }
         else
         {
            // create imposed value
            vs = double_manager->insert(keys, value_creator);
            vs->setVal(cst);
            fixed_value_creator(keys, vs);
         }

         // extra relation in the specific case of Ur imposed with r=0
         if (x_max)
         {
            assert(master.size() < 2);
            xValue<double>* other = nullptr;
            if (master.size()) other = master[0];

            keys = (col_id) ? ((col_id > 1) ? keyz : keyy) : keyx;

            if (*x_max > *y_max)
            {
               vgama = x_direction;
               cst = csteq * val;
            }
            else
            {
               vgama = y_direction;
               cst = 0.;
            }

            coeffs.resize(2);
            master.resize(2);
            int k = 0;

            // if first coefficiant of new linear relation is not null
            if (fabs(vgama(slv1)) > epsilon)
            {
               coeffs[0] = -vgama(slv1) * csteq;
               if (slv1 == slvb)
                  master[0] = vs;
               else if (other)
                  master[0] = other;
               else
                  master[0] = double_manager->insert((slv1) ? ((slv1 > 1) ? keyz : keyy) : keyx, value_creator);
               ++k;
            }
            else
               coeffs.resize(1);

            // if second coefficiant of new linear relation is not null
            if (fabs(vgama(slv2)) > epsilon)
            {
               coeffs[k] = -vgama(slv2) * csteq;
               if (slv2 == slvb)
                  master[k] = vs;
               else if (other)
                  master[k] = other;
               else
                  master[k] = double_manager->insert((slv2) ? ((slv2 > 1) ? keyz : keyy) : keyx, value_creator);
               ++k;
            }
            else
               coeffs.resize(k);

            // depending if linear relation reduce to constant or not generate last value
            if (k)
            {
               // create slave value
               vs = static_cast<xValue<double>*>(new xValueLinearCombination<double>(coeffs, master, cst));
               double_manager->insert(keys, vs);
            }
            else
            {
               // create imposed value
               vs = double_manager->insert(keys, value_creator);
               vs->setVal(cst);
               fixed_value_creator(keys, vs);
            }
         }
      }
   }
   return;
}

/*!
 *
 * Needed
 * FIELD member :
 *     xValueManagerDist<double> * getValueManager()
 *     std::vector<spacePtr>::iterator begin()
 *     std::vector<spacePtr>::iterator end()
 *
 * CREA :
 *     could be a object function, a boost or std function, a lambda function, ....
 *     In all cases must respect following signature :
 *             xValue<VT>* operator()(const xValKey& key, xValue<VT>* v)
 */
template <class FIELD, class CREA, class ITER>
void DeclareState(const FIELD& fct, const CREA& state_creator, ITER it, ITER end)
{
   auto double_manager = fct.getValueManager();
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      xFiniteElement FEM;
      FEM.setKeys(e, fct.begin(), fct.end());
      typename FIELD::value_container_t vals;
      double_manager->getValPtr(FEM.beginKey(), FEM.endKey(), vals);
      std::transform(FEM.beginKey(), FEM.endKey(), vals.begin(), vals.begin(), state_creator);
   }

   return;
}
template <class FIELD, class CREA, class ITER>
void DeclareState2(const FIELD& fct, const CREA& state_creator, ITER it, ITER end)
{
   auto double_manager = fct.getValueManager();

   std::vector<xValKey> glob;
   std::vector<xValKey> loc;
   auto space_range = xtool::make_range(fct.begin(), fct.end());
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      glob.clear();
      for (auto space : space_range)
      {
         loc.clear();
         space->getKeys(e, &loc);
         glob.insert(glob.end(), loc.begin(), loc.end());
      }
      typename FIELD::value_container_t vals;
      double_manager->getValPtr(glob.begin(), glob.end(), vals);
      std::transform(glob.begin(), glob.end(), vals.begin(), vals.begin(), state_creator);
   }

   return;
}

template <class FIELD, class ITER>
void DeleteState(const FIELD& fct, ITER it, ITER end)
{
   auto double_manager = fct.getValueManager();
   for (; it != end; ++it)
   {
      AOMD::mEntity* e = *it;
      xFiniteElement FEM;
      FEM.setKeys(e, fct.begin(), fct.end());
      typename FIELD::value_container_t vals;
      double_manager->getValPtr(FEM.beginKey(), FEM.endKey(), vals);
      std::for_each(vals.begin(), vals.end(), std::mem_fun(&xValue<typename FIELD::value_t>::delState));
   }

   return;
}

// L2Projection
template <class FIELD, class T, class ASSEMBLER, class SOLVER, class ITER>
void L2Projection(FIELD& field, const xEval<T>& given, ASSEMBLER& assembler, const xIntegrationRule& integration_rule,
                  SOLVER& solver, ITER it, ITER end)
{
   // condition to be satisfied before calling L2 projection.
   //      variables defined
   ITER beg = it;

   auto double_manager = field.getValueManager();
   xStateDofCreator<> snh(*double_manager, "l2_dofs");
   DeclareState(field, snh, beg, end);

   // This is not generic now (09/2011) as matrix container are now proposing 2
   // interfaces :
   //    * the old way : xlinalg::xCSRMatrix type where constructor use number of dof and
   //                    structure of the sparce matrix is done on the fly
   //    * the new way : xlinalg::xGenericSparseMatrix type where constructor use a
   //                    graph discribing structure of the sparce matrix
   //  here we use the old matrix container style
   typename ASSEMBLER::vector_type b(double_manager->size("l2_dofs"));
   typename ASSEMBLER::vector_type sol(double_manager->size("l2_dofs"));
   typename ASSEMBLER::matrix_type A(double_manager->size("l2_dofs"));

   assembler.setTarget(A, b);

   xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<T>>, xValOperator<xtool::xIdentity<T>>> l2;
   xFormLinearWithLoad<xValOperator<xtool::xIdentity<T>>, xEval<T>> l2_rhs(given);
   Assemble(l2, assembler, integration_rule, field, field, beg, end);
   Assemble(l2_rhs, assembler, integration_rule, field, beg, end);

   solver.connectMatrix(A);
   solver.solve(b, sol);
   // system.Solve(b, sol);
   Visit(xWriteSolutionVisitor<typename ASSEMBLER::vector_type>(sol.begin()), double_manager->begin("l2_dofs"),
         double_manager->end("l2_dofs"));

   // delete states and set values now to fixed
   DeleteState(field, beg, end);
   //  DeclareState(field, xStateFixedCreator(), beg, end);

   double_manager->clear_subset("l2_dofs");
}

// helper function for generateMeshToMeshConstrainWithL2Projection
int setDofNum(xValue<double>* v);

// generateMeshToMeshConstrainWithL2Projection
template <typename S, typename F, typename AM, typename ARHS, typename VM, typename C, typename FT, typename G, typename ITER1,
          typename ITER2>
void generateMeshToMeshConstrainWithL2Projection(const xIntegrationRule& integration_rule, S& linear_solver, const F& field,
                                                 AM& assembler_mass, ARHS& assembler_rhs, VM& val_manager, C& creator_reg,
                                                 ITER1 it_slave_elem, ITER1 it_slave_elem_end, ITER2 it_exclude_elem,
                                                 ITER2 it_exclude_elem_end, xEntityToEntity slave2master,
                                                 double thresold = 1.e-12)
{
   int numdof_master, numdof_slave, numdof;
   ITER1 it_elem_beg_slave = it_slave_elem;
   ITER2 it_elem_beg_exclude = it_exclude_elem;
   xFormBilinearWithoutLaw<xValOperator<xtool::xIdentity<FT>>, xValOperator<xtool::xIdentity<FT>>> l2_bilin;
   xIntegrateFormCommand command(&l2_bilin);
   std::map<ITER1, std::vector<xValue<typename F::value_t>*>> slave_vals;
   std::map<ITER1, xFiniteElement*> slave_FEM;
   std::map<AOMD::mEntity*, std::vector<xValue<double>*>> master_vals;
   std::map<AOMD::mEntity*, xFiniteElement*> master_FEM;
   std::vector<xValKey> slave_keys;
   typename xValueLinearCombination<typename F::value_t>::xvalues_t master_real_value;
   VM local_master_val_manager;
   VM local_slave_val_manager;

   // loop on slave to generate slave and master dofs
   for (it_slave_elem = it_elem_beg_slave; it_slave_elem != it_slave_elem_end; ++it_slave_elem)
   {
      // local
      std::vector<xValue<double>*> vals;
      xFiniteElement* FEM = new xFiniteElement;

      FEM->setKeysAndFcts((*it_slave_elem), field.begin(), field.end(), "slave");

      AOMD::mEntity* master_elem = slave2master(*it_slave_elem);

      if (!master_elem)
      {
         std::cout << "File " << __FILE__ << " line " << __LINE__ << " Master element is not attached to slave element "
                   << std::endl;
         throw;
      }
      // generate master key if not yet done
      else if (master_FEM.find(master_elem) == master_FEM.end())
      {
         // local
         std::vector<xValue<double>*> valsm;
         xFiniteElement* FEMm = new xFiniteElement;

         // finding keys for master element "master"
         FEMm->setKeysAndFcts(master_elem, field.begin(), field.end(), "master");

         // loop on keys
         xFiniteElement::iterKey itkem = FEMm->endKey("master");
         for (xFiniteElement::iterKey itkm = FEMm->beginKey("master"); itkm != itkem; ++itkm)
         {
            xValue<double>*vm, *v_true;
            xValKey& keym = (*itkm);

            // look if not already created in val manager
            // and if not, generate a std value in val manager
            if (!(v_true = static_cast<xValue<double>*>(val_manager.insert(keym, creator_reg))))
            {
               std::cout << "File " << __FILE__ << " line " << __LINE__
                         << " It appends !!! generation of a xValue<T> failed for master key in val manager" << keym << std::endl;
               throw;
            }

            // generate std value in local val manager
            if (!(vm = static_cast<xValue<double>*>(local_master_val_manager.insert(keym, creator_reg))))
            {
               std::cout << "File " << __FILE__ << " line " << __LINE__
                         << " It appends !!! generation of a xValue<T> failed for master key in local val manager " << (*itkm)
                         << std::endl;
               throw;
            }

            // set state of std value as dof using numdof
            numdof = local_master_val_manager.size("master") + 1;
            if (!vm->getState())
            {
               local_master_val_manager.add(vm, "master");
               vm->setState(new xStateOfValueDof(numdof));
               master_real_value.push_back(v_true);
            }
            // store value
            valsm.push_back(vm);
         }

         // store for loop below
         master_FEM.insert(std::make_pair(master_elem, FEMm));
         master_vals.insert(std::make_pair(master_elem, valsm));
      }

      // create slave values in ValManager
      xFiniteElement::iterKey itke = FEM->endKey("slave");
      for (xFiniteElement::iterKey itk = FEM->beginKey("slave"); itk != itke; ++itk)
      {
         xValue<double>* v;
         xValKey& key = (*itk);

         // check that this key doesn't already exist in val_manager
         if ((v = val_manager.find(key)))
         {
            // if the key is related to a entity not excluded there is a problem
            // otherwise no
            if (std::find(it_elem_beg_exclude, it_exclude_elem_end, key.getEnti()) == it_exclude_elem_end)
            {
               std::cout << "File " << __FILE__ << " line " << __LINE__ << " key " << key
                         << " already exist in val manger ?\n You can't introduce a MPC on that value" << std::endl;
               throw;
            }
         }
         // if not already created in local slave val manager insert std value
         if (!(v = static_cast<xValue<double>*>(local_slave_val_manager.insert(key, creator_reg))))
         {
            std::cout << "File " << __FILE__ << " line " << __LINE__
                      << " It appends !!! generation of a xValue<T> failed for slave key in val manager" << key << std::endl;
            throw;
         }
         numdof = local_slave_val_manager.size("slave") + 1;
         if (!v->getState())
         {
            local_slave_val_manager.add(v, "slave");
            v->setState(new xStateOfValueDof(numdof));
            slave_keys.push_back(key);
         }
         vals.push_back(v);
      }

      // store for loop below
      slave_FEM.insert(make_pair(it_slave_elem, FEM));
      slave_vals.insert(make_pair(it_slave_elem, vals));
   }

   // number of dofs
   numdof_slave = local_slave_val_manager.size("slave");
   numdof_master = local_master_val_manager.size("master");

   // setting RHS container
   // it follow xlinalg::xDenseMatrix paradigm.
   typename ARHS::matrix_type B(numdof_slave, numdof_master);
   assembler_rhs.setTarget(B);

   // setting solution container
   // it follow xlinalg::xDenseMatrix paradigm as B
   typename ARHS::matrix_type SOL(numdof_slave, numdof_master);

   // loop on slave to generate matrix graph
   G graph(numdof_slave, numdof_slave, 0);//Taucs is gone, graph is unsym (cause we use SuperLU)
   for (it_slave_elem = it_elem_beg_slave; it_slave_elem != it_slave_elem_end; ++it_slave_elem)
   {
      std::vector<xValue<double>*>& vals_slave = slave_vals.find(it_slave_elem)->second;
      std::vector<int> dof_num(vals_slave.size());
      std::transform(vals_slave.begin(), vals_slave.end(), dof_num.begin(), setDofNum);

      int i, nb_dofs = dof_num.size();

      // treating dofs if any
      if (nb_dofs > 1)
      {
         // sort numdof in ascending order
         std::vector<int>::iterator itbeg = dof_num.begin();
         std::vector<int>::iterator itend = dof_num.end();
         sort(itbeg, itend);

         // loop on pseudo colone of the elementary matrix
         for (i = 0; i < nb_dofs; ++i)
         {
            graph.addLinesSymBlock(i, nb_dofs, dof_num[i], &dof_num[0]);
         }
      }
      else
      {
         i = dof_num[0];
         graph.add(i, i);
      }
   }

   // setting nnz
   graph.countNNZ();

   // setting matrix container
   // it follow xlinalg::xGenericSparseMatrix paradigm.
   typename AM::matrix_type A(graph);

   assembler_mass.setTarget(A);

   // for types matching generate cast
   // xAssembler & assembler_rhs_cast = ( xAssembler & )assembler_rhs;
   // xAssembler & assembler_mass_cast = ( xAssembler & )assembler_mass;

   // loop on slave to generate mass matrix and RHS
   for (it_slave_elem = it_elem_beg_slave; it_slave_elem != it_slave_elem_end; ++it_slave_elem)
   {
      // finding keys and values for slave element
      xFiniteElement* FEM_slave = slave_FEM.find(it_slave_elem)->second;
      std::vector<xValue<double>*>& vals_slave = slave_vals.find(it_slave_elem)->second;

      // integration of slave for mass matrix
      xGeomElem geo_slave((*it_slave_elem));
      l2_bilin.init(&geo_slave, FEM_slave->getFcts("slave"));
      integration_rule.accept(command, (*it_slave_elem));
      l2_bilin.finalize();
      assembler_mass.assemble_matrix(vals_slave.begin(), vals_slave.end(), l2_bilin.getMatrix());

      // get master element attached to this slave. No check as it have already been tested above
      AOMD::mEntity* master_elem = slave2master(*it_slave_elem);

      // finding keys and values for master element
      xFiniteElement* FEM_master = master_FEM.find(master_elem)->second;
      std::vector<xValue<double>*>& vals_master = master_vals.find(master_elem)->second;

      // integration of slave on master for RHS matrix
      xGeomElem geo_master(master_elem);
      l2_bilin.init(&geo_slave, FEM_slave->getFcts("slave"), &geo_master, FEM_master->getFcts("master"));
      integration_rule.accept(command, (*it_slave_elem));
      l2_bilin.finalize();
      assembler_rhs.assemble_matrix(vals_slave.begin(), vals_slave.end(), vals_master.begin(), vals_master.end(),
                                    l2_bilin.getMatrix());

      // free slave memory
      delete FEM_slave;
   }

   // free master memory
   std::map<AOMD::mEntity*, xFiniteElement*>::iterator itMFe = master_FEM.end();
   for (std::map<AOMD::mEntity*, xFiniteElement*>::iterator itMF = master_FEM.begin(); itMF != itMFe; ++itMF)
      delete ((*itMF).second);

   // solve problem
   linear_solver.connectMatrix(A);
   linear_solver.solve(B, SOL);

   // loop on slave dofs to generate MPC
   for (numdof = 0; numdof < numdof_slave; ++numdof)
   {
      xValKey& key = slave_keys[numdof];
      // filter excluded entities from MPC creation
      if (std::find(it_elem_beg_exclude, it_exclude_elem_end, key.getEnti()) == it_exclude_elem_end)
      {
         typename xValueLinearCombination<typename F::value_t>::coeffs_t coeffs;
         coeffs.reserve(numdof_master);
         typename xValueLinearCombination<typename F::value_t>::xvalues_t master_real_value_nz;
         master_real_value_nz.reserve(numdof_master);
         for (int j = 0; j < numdof_master; ++j)
         {
            const double r = SOL(numdof, j);

            if (r > thresold || r < -thresold)
            {
               coeffs.push_back(r);
               master_real_value_nz.push_back(master_real_value[j]);
            }
         }
         xValue<double>* v =
             static_cast<xValue<double>*>(new xValueLinearCombination<typename F::value_t>(coeffs, master_real_value_nz));
         val_manager.insert(key, v);
      }
   }
}

// helper recursive function for DeclareInterpolationHanging
// do not use directely
template <typename S, typename F, typename AM, typename ARHS, typename VM, typename C, typename FT, typename G>
void DeclareHangingOnEntity(const xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& isHangingOn,
                            const xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>>& downGroup,
                            const xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>>& bndGroup,
                            const xIntegrationRule& integration_rule, AOMD::mEntity* nh, std::set<AOMD::mEntity*>& treated,
                            S& linear_solver, const F& field, AM& assembler_mass, ARHS& assembler_rhs, VM& val_manager,
                            C& creator_reg)
{
   // get attached vector
   std::vector<AOMD::mEntity*>& slave = const_cast<std::vector<AOMD::mEntity*>&>(*downGroup.getData(*nh));
   std::vector<AOMD::mEntity*>& exclude = const_cast<std::vector<AOMD::mEntity*>&>(*bndGroup.getData(*nh));

   // check that exclude node which are hanging are already treated

   for (auto extrem_nod : exclude)
   {
      // if node and hanging and not treated yet
      if (extrem_nod->getLevel() == 0 && isHangingOn.getData(*extrem_nod) &&
          (treated.find(extrem_nod) ==
           treated.end()))  // extrem_nod->getAttachedEntity(id_is_hanging_on) && ( treated.find(extrem_nod) == treated.end() ) )
      {
         DeclareHangingOnEntity<S, F, AM, ARHS, VM, C, FT, G>(isHangingOn, downGroup, bndGroup, integration_rule, extrem_nod,
                                                              treated, linear_solver, field, assembler_mass, assembler_rhs,
                                                              val_manager, creator_reg);
      }
   }

   // generate MPC
   // xAttached slave2master(id_is_hanging_on);
   xEntityToEntity slave2master = [&isHangingOn](AOMD::mEntity* es) {
      AOMD::mEntity* const* ppmaster = isHangingOn.getData(*es);
      if (ppmaster) return const_cast<AOMD::mEntity*>(*ppmaster);
      AOMD::mEntity* pe = nullptr;
      return pe;
   };
   generateMeshToMeshConstrainWithL2Projection<S, F, AM, ARHS, VM, C, FT, G, std::vector<AOMD::mEntity*>::iterator,
                                               std::vector<AOMD::mEntity*>::iterator>(
       integration_rule, linear_solver, field, assembler_mass, assembler_rhs, val_manager, creator_reg, slave.begin(),
       slave.end(), exclude.begin(), exclude.end(), slave2master);

   treated.insert(nh);
}

/*!
 * DeclareInterpolationHanging
 *
 * Needed
 * F member :
 *     std::vector<spacePtr>::iterator begin()
 *     std::vector<spacePtr>::iterator end()
 *
 * C :
 *     could be a object function, a boost or std function, a lambda function, ....
 *     In all cases must respect following signature :
 *         xValue<typename FIELD::value_t>* (const xValKey&) const
 */
template <typename ITER, typename S, typename F, typename AM, typename ARHS, typename VM, typename C, typename FT, typename G>
void DeclareInterpolationHanging(const int dim, const xIntegrationRule& integration_rule, ITER it, ITER itend, S& linear_solver,
                                 const F& field, AM& assembler_mass, ARHS& assembler_rhs, VM& val_manager, C& creator_reg,
                                 const xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& isHangingOn,
                                 const xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>>& downGroup,
                                 const xinterface::aomd::xAttachedDataManagerAOMD<std::vector<AOMD::mEntity*>>& bndGroup)
{
   std::set<AOMD::mEntity*> treated;

   // loop on node to treat all hanging entities
   for (; it != itend; ++it)
   {
      // set as a pointer
      AOMD::mEntity* nh = (*it);

      // if hanging
      // if (nh->getAttachedEntity(id_is_hanging_on))
      if (isHangingOn.getData(*nh))
      {
         // if this entity is not treated yet, do so
         if (treated.find(nh) == treated.end())
         {
            // DeclareHangingOnEntity < S,F,AM,ARHS,VM,C,FT,G >(id_is_hanging_on,id_down_groupe,id_bnd_groupe,integration_rule,
            // nh,
            //                                                 treated,linear_solver, field, assembler_mass, assembler_rhs,
            //                                                 val_manager, creator_reg);

            DeclareHangingOnEntity<S, F, AM, ARHS, VM, C, FT, G>(isHangingOn, downGroup, bndGroup, integration_rule, nh, treated,
                                                                 linear_solver, field, assembler_mass, assembler_rhs, val_manager,
                                                                 creator_reg);
         }
      }
   }
}

/// fills a xFieldPointwise with values coming from an evaluator
/// can be used to init a field, or copy one into another
/// The name FillField (and not FillFieldpointwise) reflects the fact that
/// one may code a version for approximations (xField)  using an L2 projection

// void FillField(xFieldPointwise& fct,const xEval<T> &given, const xIntegrationRule& integration_rule,ITER it, ITER end)

template <typename T, class ITER>
void FillField(xFieldPointwise& fct, const xEval<T>& given, const xIntegrationRule& integration_rule, ITER it, ITER end,
               xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   xSetFieldPointwiseVisitor<T> Init(fct.valstr(), given);
   xPointwiseCommand init_command(Init, **(fct.begin()), *fct.getValueManager());
   ApplyCommandOnIntegrationRule(init_command, integration_rule, it, end, integ2appro);
}

template <typename T, class ITER>
double FillFieldAndCheckChanges(xFieldPointwise& fct, const xEval<T>& given, const xIntegrationRule& integration_rule, ITER it,
                                ITER end, xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   xSetFieldPointwiseCheckChangeVisitor<T> Init(fct.valstr(), given);
   xPointwiseCommand init_command(Init, **(fct.begin()), *fct.getValueManager());
   ApplyCommandOnIntegrationRule(init_command, integration_rule, it, end, integ2appro);
   return Init.getChangeRatio();
   // return Init.getNumberOfChanges();
}

/// Add values to already existing field
template <typename T, class ITER>
void AddToField(xFieldPointwise& fct, const xEval<T>& given, const xIntegrationRule& integration_rule, ITER it, ITER end,
                xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   xAddToFieldPointwiseVisitor<T> Init(fct.valstr(), given);
   xPointwiseCommand init_command(Init, **(fct.begin()), *fct.getValueManager());
   ApplyCommandOnIntegrationRule(init_command, integration_rule, it, end, integ2appro);
}

template <class FIELD, class ITER>
int ComputeNnz(FIELD& field, ITER start, ITER last)
{
   size_t dofsize = field.getValueManager()->size("dofs");
   vector<std::set<size_t>> graph(dofsize);
   // -- 1. Graph construction -------------------------------------------------
   for (; start != last; ++start)
   {
      xFiniteElement fem;
      fem.setKeys(*start, field.begin(), field.end());

      vector<xValKey>* keys = fem.getKeys();

      for (vector<xValKey>::const_iterator itk1 = keys->begin(); itk1 != keys->end(); ++itk1)
      {
         xValue<double>* v1 = field.getValueManager()->find(*itk1);
         vector<xValKey>::const_iterator itk2 = itk1;

         // if( itk2 != keys->end() ) ++itk2;

         for (; itk2 != keys->end(); ++itk2)
         {
            xValue<double>* v2 = field.getValueManager()->find(*itk2);

            if (v1 != nullptr && v2 != nullptr)
            {
               xStateOfValueDof* s1 = dynamic_cast<xStateOfValueDof*>(v1->getState());
               xStateOfValueDof* s2 = dynamic_cast<xStateOfValueDof*>(v2->getState());

               if (s1 != nullptr && s2 != nullptr)
               {
                  graph[s1->Numdof - 1].insert(s2->Numdof - 1);
                  graph[s2->Numdof - 1].insert(s1->Numdof - 1);
               }
            }
         }
      }
   }
   vector<std::set<size_t>>::const_iterator it = graph.begin();
   vector<std::set<size_t>>::const_iterator ite = graph.end();
   int nnz = 0;
   for (; it != ite; ++it)
   {
      nnz += it->size();
   }
   return nnz;
}

template <class FIELD, class ITER>
int ComputeGraphAndNnz(FIELD& field, ITER start, ITER last, vector<std::set<size_t>>& graph)
{
   size_t dofsize = field.getValueManager()->size("dofs");
   graph.resize(dofsize);
   // -- 1. Graph construction -------------------------------------------------
   for (; start != last; ++start)
   {
      xFiniteElement fem;
      fem.setKeys(*start, field.begin(), field.end());

      vector<xValKey>* keys = fem.getKeys();

      for (vector<xValKey>::const_iterator itk1 = keys->begin(); itk1 != keys->end(); ++itk1)
      {
         xValue<double>* v1 = field.getValueManager()->find(*itk1);
         vector<xValKey>::const_iterator itk2 = itk1;

         // if( itk2 != keys->end() ) ++itk2;

         for (; itk2 != keys->end(); ++itk2)
         {
            xValue<double>* v2 = field.getValueManager()->find(*itk2);

            if (v1 != nullptr && v2 != nullptr)
            {
               xStateOfValueDof* s1 = dynamic_cast<xStateOfValueDof*>(v1->getState());
               xStateOfValueDof* s2 = dynamic_cast<xStateOfValueDof*>(v2->getState());

               if (s1 != nullptr && s2 != nullptr)
               {
                  graph[s1->Numdof - 1].insert(s2->Numdof - 1);
                  graph[s2->Numdof - 1].insert(s1->Numdof - 1);
               }
            }
         }
      }
   }
   vector<std::set<size_t>>::const_iterator it = graph.begin();
   vector<std::set<size_t>>::const_iterator ite = graph.end();
   int nnz = 0;
   for (; it != ite; ++it)
   {
      nnz += it->size();
   }
   return nnz;
}

template <class FIELD, class ITER>
void TreatmentOfEssentialEnv(const FIELD& fct, ITER itenv, ITER endenv, const std::vector<string>& vars, xMesh* mesh)
{
   //	const bool debug = xdebug_flag;
   for (; itenv != endenv; ++itenv)
   {
      const xEnv& env = *itenv;
      string phys = env.Phys;
      int type = env.Type;
      int entity = env.Entity;
      if (find(vars.begin(), vars.end(), phys) != vars.end())
      {
         if (type == FIX)
         {
            xClassRegion bc(mesh, entity, env.getDimension());
            double val = env.getValue();
            DirichletBoundaryCondition(fct, phys, bc.begin(), bc.end(), val);
            //				DirichletBoundaryCondition (fct, xIsPhys(phys), bc.begin(), bc.end(), val);
         }
      }
   }
   return;
}

template <class FIELD>
void TreatmentOfEssentialEnv(const FIELD& fct, const std::vector<string>& vars, xMesh* mesh,
                             std::list<xEntityFilter> dirichlet_filters)
{
   xRegion all(mesh);
   xAcceptUnion filter_dirichlet(dirichlet_filters);
   xFilteredRegion<xIter, xAcceptUnion> dir_reg(all.begin(1), all.end(1), filter_dirichlet);

   // const bool debug = xdebug_flag;
   for (std::vector<string>::const_iterator it = vars.begin(); it != vars.end(); it++)
   {  // iterateur sur vector strings
      string phys = (*it);
      double val = xtool::xDataType<typename FIELD::value_t>::zero();
      DirichletBoundaryCondition(fct, phys, dir_reg.begin(), dir_reg.end(), val);
   }
   return;
}

template <class FIELD>
void TreatmentOfEssentialEnvIntersection(const FIELD& fct, const std::vector<string>& vars, xMesh* mesh,
                                         std::list<xEntityFilter> dirichlet_filters)
{
   xRegion all(mesh);
   xAcceptIntersection filter_dirichlet(dirichlet_filters);
   xFilteredRegion<xIter, xAcceptIntersection> dir_reg(all.begin(1), all.end(1), filter_dirichlet);

   // const bool debug = xdebug_flag;
   for (std::vector<string>::const_iterator it = vars.begin(); it != vars.end(); it++)
   {  // iterateur sur vector strings
      string phys = (*it);
      double val = xtool::xDataType<typename FIELD::value_t>::zero();
      DirichletBoundaryCondition(fct, phys, dir_reg.begin(), dir_reg.end(), val);
   }
   return;
}

/*!
 *
 * Needed
 * see DeclareInterpolationRigid
 */
template <class FIELD, class CREA>
void DeclareInterpolationRigid(const FIELD& fct, const CREA& value_creator, const std::vector<string>& vars, xMesh* mesh,
                               std::list<xEntityFilter> dirichlet_filters)
{
   xRegion all(mesh);
   xAcceptUnion filter_dirichlet(dirichlet_filters);
   xFilteredRegion<xIter, xAcceptUnion> dir_reg(all.begin(1), all.end(1), filter_dirichlet);
   for (std::vector<string>::const_iterator it = vars.begin(); it != vars.end(); it++)
   {  // iterateur sur vector strings
      string phys = (*it);
      DeclareInterpolationRigid(fct, value_creator, phys, dir_reg.begin(), dir_reg.end());
   }
   return;
}

template <class ASSEMBLER, class FIELD, class ITER>
void AssembleNaturalEnvScalar(ASSEMBLER& assembler, const xIntegrationRule& integration_rule, const FIELD& test, ITER itenv,
                              ITER endenv, const std::vector<string>& vars, xMesh* mesh,
                              xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   const bool debug = xdebug_flag;
   for (; itenv != endenv; ++itenv)
   {
      const xEnv& env = *itenv;
      std::vector<string>::const_iterator it = find(vars.begin(), vars.end(), env.Phys);
      if (it != vars.end())
      {
         assert(env.Type == FIX);
         double val = env.getValue();
         if (debug) std::cout << "val is " << val << std::endl;
         xEvalConstant<typename FIELD::value_t> flux(val);
         xFormLinearWithLoad<xValOperator<xtool::xIdentity<double>>, xEvalConstant<typename FIELD::value_t>,
                             typename FIELD::value_t>
             lin(flux);  // Shape functions are assumed double !
         xClassRegion bc(mesh, env.Entity, env.getDimension());
         Assemble(lin, assembler, integration_rule, test, bc.begin(), bc.end(), integ2appro);
      }
   }
   return;
}

template <class ASSEMBLER, class FIELD, class ITER>
void AssembleNaturalEnvVector(ASSEMBLER& assembler, const xIntegrationRule& integration_rule, const FIELD& test, ITER itenv,
                              ITER endenv, const std::vector<string>& vars, xMesh* mesh,
                              xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   const bool debug = xdebug_flag;
   using xvector_type = typename xtensor::xVector<typename FIELD::value_t>;
   for (; itenv != endenv; ++itenv)
   {
      const xEnv& env = *itenv;
      std::vector<string>::const_iterator it = find(vars.begin(), vars.end(), env.Phys);
      if (it != vars.end())
      {
         assert(env.Type == FIX);
         xtensor::xVector<> val;
         int c = distance(vars.begin(), it);  // Warning :  this is dangerous !  it suppose that vars can only
         // have 3 string, one for each component of the vector to fix !
         val(c) = env.getValue();
         if (debug) std::cout << "val is " << val(0) << " " << val(1) << " " << val(2) << std::endl;
         xEvalConstant<xvector_type> flux(val);
         //            xFormLinearWithLoad<xValOperator<xtool::xIdentity<xtensor::xVector<> > >, xEvalConstant<xvector_type>,
         //            typename FIELD::value_t > lin(flux);//Shape functions are assumed double, the assembled quantity is of type
         //            FIELD::value_t
         xFormLinearWithLoad<xValOperator<xtool::xIdentity<xtensor::xVector<>>>, xEvalConstant<xvector_type>> lin(
             flux);  // Shape functions are assumed double, the assembled quantity is of type FIELD::value_t
         xClassRegion bc(mesh, env.Entity, env.getDimension());
         Assemble(lin, assembler, integration_rule, test, bc.begin(), bc.end(), integ2appro);
      }
   }
   return;
}

string GetFileNameWithProcessRank(const string& base, const string& exte);

template <typename VT>
void addLinCombDOF(const xStateOfValueLinearCombination<VT>* equ, std::vector<int>& dof_num, std::vector<int>::size_type& m)
{
   xStateOfValue* s;
   xStateOfValue::state_type state;

   // number of equation termes
   size_t nb_equ = equ->size();

   // reseve whatever append
   if (dof_num.capacity() - m < nb_equ)
   {
      m += nb_equ;
      dof_num.reserve(m);
   }

   // loop on equation termes
   for (size_t k = 0; k < nb_equ; ++k)
   {
      s = equ->state(k);
      state = s->state;
      if (state == xStateOfValue::DOF)
         dof_num.push_back((static_cast<const xStateOfValueDof*>(s))->Numdof - 1);
      else if (state == xStateOfValue::LINEARCOMBINATION)
         addLinCombDOF((static_cast<const xStateOfValueLinearCombination<VT>*>(s)), dof_num, m);
   }
}

// This methode, giving a field and elements generate graph structure
// Here symetrique pattern is considered as trial/test are identical
template <typename G, typename FIELD, typename ITER>
void AssembleGraph(G& graph, FIELD& field, ITER start, ITER last,
                   xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   // local
   int i, nb_dofs;
   std::vector<int>::size_type l;
   xValue<typename FIELD::value_t>* v;

   // get double manager associated to this field
   auto dbmng = field.getValueManager();

   // loop on entity
   for (; start != last; ++start)
   {
      // get keys for curent entity
      AOMD::mEntity* e_appro = integ2appro(*start);
      xFiniteElementKeysOnly fem;
      fem.setKeys(e_appro, field.begin(), field.end());

      std::vector<xValKey>* keys = fem.getKeys();
      std::vector<xValKey>::const_iterator itkend = keys->end();
      std::vector<int>::size_type nb_keys = keys->size();

      // create container which will store number of identified DOF directly or from LinearCombination
      std::vector<int> dof_num;

      // pres allocating container with  number of possible ddls from keys
      dof_num.reserve(nb_keys);

      // loop on keys to identify DOF and LinearCombination ddls which are the only ones which give
      // non zero term in the matrix
      for (std::vector<xValKey>::const_iterator itk = keys->begin(); itk != itkend; ++itk)
      {
         // see if value exist in double manager
         v = dbmng->find(*itk);

         // if yes store
         if (v != nullptr)
         {
            // get state
            xStateOfValue* s = v->getState();
            xStateOfValue::state_type state = s->state;
            l = dof_num.size();
            // check state
            if (state == xStateOfValue::DOF)
               dof_num.push_back((static_cast<const xStateOfValueDof*>(s))->Numdof - 1);
            else if (state == xStateOfValue::LINEARCOMBINATION)
               addLinCombDOF((static_cast<const xStateOfValueLinearCombination<typename FIELD::value_t>*>(s)), dof_num, l);
         }
      }

      // treating dofs if any
      nb_dofs = dof_num.size();
      if (nb_dofs > 1)
      {
         // sort numdof in ascending order
         std::vector<int>::iterator itbeg = dof_num.begin();
         std::vector<int>::iterator itend = dof_num.end();
         sort(itbeg, itend);

         // eliminate duplicate entry if any : during recursive process of addLinCombDOF somme dof may appears many times
         // i.e. a dof may be free and be a terme of one linear cobination and be a terme of a other linear cobination induced by
         // the first one and so one
         std::vector<int>::iterator itendu = std::unique(itbeg, itend);
         nb_dofs = itendu - itbeg;

         if (nb_dofs > 1)
         {
            // loop on pseudo colone of the elementary matrix
            for (i = 0; i < nb_dofs; ++i)
            {
               graph.addLinesSymBlock(i, nb_dofs, dof_num[i], &dof_num[0]);
            }
         }
         else
         {
            i = dof_num[0];
            graph.add(i, i);
         }
      }
      else if (nb_dofs)
      {
         i = dof_num[0];
         graph.add(i, i);
      }
   }

   // resetting nnz as term may have been added here
   graph.countNNZ();
}

// This methode, giving 2 fields and elements generate graph structure
// Here unsymetrique pattern is considered as trial/test are different
template <typename G, typename FIELD, typename ITER>
void AssembleGraph(G& graph, FIELD& test, FIELD& trial, ITER start, ITER last,
                   xEntityToEntity integ2appro_test = xtool::xIdentity<AOMD::mEntity*>(),
                   xEntityToEntity integ2appro_trial = xtool::xIdentity<AOMD::mEntity*>())
{
   // eliminate the case where user doesn't pay attention to the fact that it existe a methode
   // when test=trial
   // if (&test == &trial && std::ref(integ2appro_test) == std::ref(integ2appro_trial) )
   if (&test == &trial)
   {
      AssembleGraph(graph, test, start, last, integ2appro_test);
      return;
   }

   // local
   int i, h, nb_dofs_test, nb_dofs_trial;
   std::vector<int>::size_type l;
   xValue<typename FIELD::value_t>* v;
   bool sym = graph.isSym();

   // get double manager associated to  fields
   auto dbmng_test = test.getValueManager();
   auto dbmng_trial = trial.getValueManager();

   // loop on entity
   for (; start != last; ++start)
   {
      // create container which will store number of identified DOF directly or from LinearCombination
      std::vector<int> dof_num_test;
      std::vector<int> dof_num_trial;

      // get keys for curent entity
      AOMD::mEntity* e_appro_test = integ2appro_test(*start);
      AOMD::mEntity* e_appro_trial = integ2appro_trial(*start);
      xFiniteElementKeysOnly fem;
      fem.setKeys(e_appro_test, test.begin(), test.end(), "test");

      // create stuff for test
      std::vector<xValKey>* keys = fem.getKeys("test");
      std::vector<xValKey>::const_iterator itkend = keys->end();
      std::vector<int>::size_type nb_keys_test = keys->size();

      // pres allocating container with  number of possible ddls from keys
      dof_num_test.reserve(nb_keys_test);

      // loop on keys to identify DOF and LinearCombination ddls which are the only ones which give
      // non zero term in the matrix
      for (std::vector<xValKey>::const_iterator itk = keys->begin(); itk != itkend; ++itk)
      {
         // see if value exist in double manager
         v = dbmng_test->find(*itk);

         // if yes store
         if (v != nullptr)
         {
            xStateOfValue* s = v->getState();
            xStateOfValue::state_type state = s->state;
            l = dof_num_test.size();
            if (state == xStateOfValue::DOF)
               dof_num_test.push_back((static_cast<const xStateOfValueDof*>(s))->Numdof - 1);
            else if (state == xStateOfValue::LINEARCOMBINATION)
               addLinCombDOF((static_cast<const xStateOfValueLinearCombination<typename FIELD::value_t>*>(s)), dof_num_test, l);
         }
      }

      // sort dof_num_test in ascending order
      std::vector<int>::iterator itbeg = dof_num_test.begin();
      std::vector<int>::iterator itend = dof_num_test.end();
      std::sort(itbeg, itend);

      // eliminate duplicate entry if any : during recursive process of addLinCombDOF somme dof may appears many times
      // i.e. a dof may be free and be a terme of one linear cobination and be a terme of a other linear cobination induced by the
      // first one and so one
      std::vector<int>::iterator itendu = std::unique(itbeg, itend);
      nb_dofs_test = itendu - itbeg;

      // create stuff for trial
      fem.setKeys(e_appro_trial, trial.begin(), trial.end(), "trial");
      keys = fem.getKeys("trial");
      itkend = keys->end();
      std::vector<int>::size_type nb_keys_trial = keys->size();

      // pres allocating container with  number of possible ddls from keys
      dof_num_trial.reserve(nb_keys_trial);

      // loop on keys to identify DOF and LinearCombination ddls which are the only ones which give
      // non zero term in the matrix
      for (std::vector<xValKey>::const_iterator itk = keys->begin(); itk != itkend; ++itk)
      {
         // see if value exist in double manager
         v = dbmng_trial->find(*itk);

         // if yes store
         if (v != nullptr)
         {
            xStateOfValue* s = v->getState();
            xStateOfValue::state_type state = s->state;
            l = dof_num_trial.size();
            if (state == xStateOfValue::DOF)
               dof_num_trial.push_back((static_cast<const xStateOfValueDof*>(s))->Numdof - 1);
            else if (state == xStateOfValue::LINEARCOMBINATION)
               addLinCombDOF((static_cast<const xStateOfValueLinearCombination<typename FIELD::value_t>*>(s)), dof_num_trial, l);
         }
      }

      // sort dof_num_trial in ascending order
      itbeg = dof_num_trial.begin();
      itend = dof_num_trial.end();
      std::sort(itbeg, itend);

      // eliminate duplicate entry if any : during recursive process of addLinCombDOF somme dof may appears many times
      // i.e. a dof may be free and be a terme of one linear cobination and be a terme of a other linear cobination induced by the
      // first one and so one
      itendu = std::unique(itbeg, itend);
      nb_dofs_trial = itendu - itbeg;

      // treating dofs if any
      if (nb_dofs_test > 1 && nb_dofs_trial > 1)
      {
         // check
         // for sym=true this block must be in the lower part of the global
         // matrix => dof number of test must be higher then trial ones
         // if not they have to be swaped some how to get transposed termes
         //
         // if incorrect numbering i.e at least one dof number of test is lower
         // then trial ones
         if (sym && (dof_num_trial[nb_dofs_trial - 1] > dof_num_test[0]))
         {
            // simplest case : all terms have to be swaped
            if (dof_num_trial[0] > dof_num_test[nb_dofs_test - 1])
            {
               // loop on pseudo colones of the elementary matrix with test insted of trial
               for (i = 0; i < nb_dofs_test; ++i)
               {
                  graph.addLinesUnsymBlock(nb_dofs_trial, dof_num_test[i], &dof_num_trial[0]);
               }
            }
            // only partial swap
            // use of  special loop
            else
            {
               std::vector<int>::iterator it_beg = dof_num_test.begin();
               std::vector<int>::iterator it_end = dof_num_test.end();
               std::vector<int>::iterator itlower;

               // loop on trial with upper test
               for (i = 0; i < nb_dofs_trial; ++i)
               {
                  itlower = lower_bound(it_beg, it_end, dof_num_trial[i]);
                  h = it_end - itlower;
                  if (h) graph.addLinesUnsymBlock(h, dof_num_trial[i], &(*itlower));
               }

               it_beg = dof_num_trial.begin();
               it_end = dof_num_trial.end();

               // loop on test with upper trial
               for (i = 0; i < nb_dofs_test; ++i)
               {
                  itlower = lower_bound(it_beg, it_end, dof_num_test[i]);
                  h = it_end - itlower;
                  if (h) graph.addLinesUnsymBlock(h, dof_num_test[i], &(*itlower));
               }
            }
         }
         // if not sym or sym with correct numbering
         else
         {
            // loop on pseudo colones of the elementary matrix
            for (i = 0; i < nb_dofs_trial; ++i)
            {
               graph.addLinesUnsymBlock(nb_dofs_test, dof_num_trial[i], &dof_num_test[0]);
            }
         }
      }
      else if (nb_dofs_test && nb_dofs_trial > 1)
      {
         h = dof_num_test[0];

         // check
         // if incorrect numbering i.e  the  dof number of test is lower
         // then at least one trial ones
         if (sym && (dof_num_trial[nb_dofs_trial - 1] > h))
         {
            // simplest case : role of test and trial have to be swaped
            if (dof_num_trial[0] > h)
            {
               graph.addLinesUnsymBlock(nb_dofs_trial, h, &dof_num_trial[0]);
            }
            else
            {
               //  test with upper trial
               std::vector<int>::iterator it_beg = dof_num_trial.begin();
               std::vector<int>::iterator it_end = dof_num_trial.end();
               std::vector<int>::iterator itlower = lower_bound(it_beg, it_end, h);
               // here it_end-itlower >0 otherwise we won't be there
               graph.addLinesUnsymBlock(it_end - itlower, h, &(*itlower));

               // loop on trial with upper test
               for (; it_beg != itlower; ++it_beg) graph.add(h, (*it_beg));
            }
         }
         else
         {
            for (i = 0; i < nb_dofs_trial; ++i) graph.add(h, dof_num_trial[i]);
         }
      }
      else if (nb_dofs_test > 1 && nb_dofs_trial)
      {
         h = dof_num_trial[0];

         // sort numdof in ascending order
         sort(dof_num_test.begin(), dof_num_test.end());
         // check
         // if incorrect numbering i.e  the  dof number of trial is greater
         // then at least one test ones
         if (sym && (h > dof_num_test[0]))
         {
            // simplest case : role of test and trial have to be swaped
            if (h > dof_num_test[nb_dofs_test - 1])
            {
               for (i = 0; i < nb_dofs_test; ++i) graph.add(h, dof_num_test[i]);
            }
            else
            {
               //  trial with upper test
               std::vector<int>::iterator it_beg = dof_num_test.begin();
               std::vector<int>::iterator it_end = dof_num_test.end();
               std::vector<int>::iterator itlower = lower_bound(it_beg, it_end, h);
               // here it_end-itlower >0 otherwise we won't be there
               graph.addLinesUnsymBlock(it_end - itlower, h, &(*itlower));

               // loop on test with upper trial
               for (; it_beg != itlower; ++it_beg) graph.add(h, (*it_beg));
            }
         }
         else
         {
            graph.addLinesUnsymBlock(nb_dofs_test, h, &dof_num_test[0]);
         }
      }
      else if (nb_dofs_test && nb_dofs_trial)
      {
         // for security as it never been tested .....
         std::cout << "In file " << __FILE__ << " line " << __LINE__ << " compiled " << __DATE__ << " at " << __TIME__
                   << std::endl;
         std::cout << "Sorry but you are the first to try this part of the methode.\n It have not been tested befor and may be "
                      "completely bugged.\n Please test,debug,validate and remove this message and the throw if ok"
                   << std::endl;
         throw;

         if (sym && dof_num_test[0] < dof_num_trial[0])
            graph.add(dof_num_trial[0], dof_num_test[0]);
         else
            graph.add(dof_num_test[0], dof_num_trial[0]);
      }
   }

   // resetting nnz as term may have been added here
   graph.countNNZ();
}

template <typename EVALSTRESS, typename IR, typename ITER>
void ComputResultant(xtensor::xVector<>& res, ITER first, ITER last, EVALSTRESS& eval_stress, IR& integ_rule,
                     xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   // create eval comand object representing sig.n
   xEvalNormal eval_normal;
   xEvalBinary<xtool::xMult<xtensor::xTensor2<>, xtensor::xVector<>, xtensor::xVector<>>> flux(eval_stress, eval_normal);
   xIntegrateEvalCommand<xEvalBinary<xtool::xMult<xtensor::xTensor2<>, xtensor::xVector<>, xtensor::xVector<>>>> integ_command(
       flux, res);

   // integrated on set of element given by user
   ApplyCommandOnIntegrationRule(integ_command, integ_rule, first, last, integ2appro);
}

template <typename EVALSTRESS, typename IR, typename ITER>
void ComputMoment(xtensor::xVector<>& res, xtensor::xPoint& orig, ITER first, ITER last, EVALSTRESS& eval_stress, IR& integ_rule,
                  xEntityToEntity integ2appro = xtool::xIdentity<AOMD::mEntity*>())
{
   // create eval comand object representing OM ^ sig.n
   xEvalNormal eval_normal;
   xEvalVectFromPoint eval_OM(orig);
   xEvalBinary<xtool::xMult<xtensor::xTensor2<>, xtensor::xVector<>, xtensor::xVector<>>> flux(eval_stress, eval_normal);
   xEvalBinary<xtool::xMod<xtensor::xVector<>, xtensor::xVector<>, xtensor::xVector<>>> moment(eval_OM, flux);
   xIntegrateEvalCommand<xEvalBinary<xtool::xMod<xtensor::xVector<>, xtensor::xVector<>, xtensor::xVector<>>>> integ_command(
       moment, res);

   // integrated on set of element given by user
   ApplyCommandOnIntegrationRule(integ_command, integ_rule, first, last, integ2appro);
}

}  // namespace xfem

#endif
