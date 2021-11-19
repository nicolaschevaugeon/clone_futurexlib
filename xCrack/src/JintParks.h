/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef _GetJint_H
#define _GetJint_H

// std
#include <iostream>
#include <map>
#include <string>
#include <vector>

// xtensor
#include "xTensor2.h"
#include "xTensor4.h"
#include "xVector.h"
// xfem
#include "xAlgorithm.h"
#include "xAssembler.h"
#include "xEval.h"
#include "xField.h"
#include "xGeomElem.h"
#include "xLevelSet.h"
#include "xRegion.h"
#include "xVariabManager.h"
#include "xVectorField.h"
// xexport
#include "xExportAlgorithm.h"
#include "xExportGmsh.h"
// xcrack
#include "CrackPostpro.h"

class lCrack;

using namespace AOMD;
// using namespace xlinalg;

namespace xcrack
{
class ValueCreatorFront
{
   typedef std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> hash_map_t1;
   typedef std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey> hash_map_t2;

  public:
   ValueCreatorFront(xfem::xValueManagerDist<double>* dd, const hash_map_t1& a, const hash_map_t1& b, const hash_map_t2& c)
       : double_manager(dd), e0d_s_loc(a), e1d_dof_length(b), e0d_e1d_dof(c)
   {
   }
   xfem::xValue<double>* operator()(const xfem::xValKey& key) const;

  private:
   xfem::xValueManagerDist<double>* double_manager;
   const hash_map_t1& e0d_s_loc;
   const hash_map_t1& e1d_dof_length;
   const hash_map_t2& e0d_e1d_dof;
};

class ValueCreatorJ3dOld
{
   typedef std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> hash_map_t1;
   typedef std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey> hash_map_t2;

  public:
   ValueCreatorJ3dOld(xfem::xValueManagerDist<double>* dd, const hash_map_t1& a, const hash_map_t1& b, const hash_map_t2& c,
                      const hash_map_t2& d)
       : double_manager(dd),
         e0d3d_s_loc(a),
         e1d_dof_length(b),
         e0d3d_e1d_dof(c),
         e0d3d_e1d_dof_other(d),
         value_creator_front(dd, e0d3d_s_loc, e1d_dof_length, e0d3d_e1d_dof),
         geom_id(xfem::xKeyInfo::getGeomId("HIERARCHICAL_1")),
         phys_id(xfem::xKeyInfo::getPhysId("J_front"))
   {
   }
   xfem::xValue<double>* operator()(const xfem::xValKey& key) const;

  private:
   xfem::xValueManagerDist<double>* double_manager;
   const hash_map_t1& e0d3d_s_loc;
   const hash_map_t1& e1d_dof_length;
   const hash_map_t2& e0d3d_e1d_dof;
   const hash_map_t2& e0d3d_e1d_dof_other;
   ValueCreatorFront value_creator_front;
   xfem::xValKey::ids_size_t geom_id;
   xfem::xValKey::ids_size_t phys_id;
};

class ValueCreatorJ3d
{
   typedef std::unordered_map<mEntity*, double, EntityHashKey, EntityEqualKey> hash_map_t1;
   typedef std::unordered_map<mEntity*, mEntity*, EntityHashKey, EntityEqualKey> hash_map_t2;

  public:
   ValueCreatorJ3d(xfem::xValueManagerDist<double>* dd, const hash_map_t1& a, const hash_map_t1& b, const hash_map_t2& c,
                   xfem::xMesh* m)
       : double_manager(dd),
         e0d3d_s_loc(a),
         e1d_dof_length(b),
         e0d3d_e1d_dof(c),
         mesh_front_dofs(m),
         value_creator_front(dd, e0d3d_s_loc, e1d_dof_length, e0d3d_e1d_dof),
         geom_id(xfem::xKeyInfo::getGeomId("HIERARCHICAL_1")),
         phys_id(xfem::xKeyInfo::getPhysId("J_front"))
   {
   }
   xfem::xValue<double>* operator()(const xfem::xValKey& key) const;

  private:
   xfem::xValueManagerDist<double>* double_manager;
   const hash_map_t1& e0d3d_s_loc;
   const hash_map_t1& e1d_dof_length;
   const hash_map_t2& e0d3d_e1d_dof;
   xfem::xMesh* mesh_front_dofs;
   ValueCreatorFront value_creator_front;
   xfem::xValKey::ids_size_t geom_id;
   xfem::xValKey::ids_size_t phys_id;
};

class xcSetParksRadialFunction : public xfem::xLevelSetModifier
{
  public:
   xcSetParksRadialFunction(const lCrack& c, const double& ro) : lsn(*c.getFieldn()), lst(*c.getFieldt()), rho(ro) {}
   void visit(xfem::xLevelSet& f, xfem::xRegion target) override;

  private:
   const xfem::xLevelSet& lsn;
   const xfem::xLevelSet& lst;
   const double rho;
};

class xcEvalFrParks : public xfem::xEval<double>
{
  public:
   xcEvalFrParks(const lCrack& c, const double& ro, const xfem::xRegion& d, const xfem::xRegion& dol)
       : crack(c), rho(ro), domain_for_integral(d), domain_for_integral_outer_layer(dol), fr_nodal_field(d)
   {
      const bool debug = true;
      xcSetParksRadialFunction set_fr(crack, rho);
      fr_nodal_field.accept(set_fr);
      if (debug)
      {
         xexport::xExportGmshAscii pexport;
         xexport::Export(fr_nodal_field, pexport, "fr_nodal_field");
      }
   }

   void operator()(const xfem::xGeomElem* geom_appro, const xfem::xGeomElem* geom_integ, double& result) const override
   {
      const bool debug = false;
      mEntity* e_appro = geom_appro->getEntity();
      xtensor::xPoint uvw_appro = geom_appro->getUVW();
      if (debug) std::cout << " xyz for fr parks " << geom_integ->getXYZ() << std::endl;
         // smooth case
#if 0
      if (!domain_for_integral_outer_layer.IsInRegion(e_appro))
	{
	  double r, theta;
	  crack.getLocalCoords(e_appro, uvw_appro, r, theta);
      result = exp(-(r*r/((r-rho)*(r-rho))));
      if (debug) std::cout << " fr smooth case " << result << std::endl;
	  return;
	}
#endif
      // interpolated case
      result = fr_nodal_field.getVal(e_appro, uvw_appro);
      if (debug) std::cout << " fr interpolated case " << result << std::endl;
      return;
   }

  private:
   const lCrack& crack;
   const double& rho;
   const xfem::xRegion& domain_for_integral;
   const xfem::xRegion& domain_for_integral_outer_layer;
   xfem::xLevelSet fr_nodal_field;
};

class xcEvalGradFrParks : public xfem::xEval<xtensor::xVector<>>
{
  public:
   xcEvalGradFrParks(const lCrack& c, const double& ro, const xfem::xRegion& d, const xfem::xRegion& dol)
       : crack(c),
         lsn(*c.getFieldn()),
         lst(*c.getFieldt()),
         rho(ro),
         domain_for_integral(d),
         domain_for_integral_outer_layer(dol),
         fr_nodal_field(d)
   {
      const bool debug = false;
      xcSetParksRadialFunction set_fr(crack, rho);
      fr_nodal_field.accept(set_fr);
      if (debug)
      {
         xexport::xExportGmshAscii pexport;
         xexport::Export(fr_nodal_field, pexport, "fr_nodal_field");
      }
   }

   void operator()(const xfem::xGeomElem* geom_appro, const xfem::xGeomElem* geom_integ,
                   xtensor::xVector<>& result) const override
   {
      const bool debug = false;
      mEntity* e_appro = geom_appro->getEntity();
      if (debug) std::cout << " xyz for fr parks " << geom_integ->getXYZ() << std::endl;
         // smooth case
#if 0
      xtensor::xPoint uvw_appro = geom_appro->getUVW();
      if (!domain_for_integral_outer_layer.IsInRegion(e_appro))
	{
	  double r, theta;
	  crack.getLocalCoords(e_appro, uvw_appro, r, theta);
	  xtensor::xVector<> e1 = lst.getGrad(e_appro, uvw_appro);
	  xtensor::xVector<> e2 = lsn.getGrad(e_appro, uvw_appro);

	  double fr = exp(-(r*r/((r-rho)*(r-rho))));
	  double grad_fr = 2*r*rho/((r-rho)*(r-rho)*(r-rho))*fr;

	  xtensor::xVector<> nr = e1*cos(theta) + e2*sin(theta); 


          if (debug) 
	    {
              std::cout << " in grad fr parks " << std::endl;
              std::cout << " r " << r << " theta " << theta  << " rho " << rho << std::endl;
              std::cout << " grad fr sca " << grad_fr << std::endl;
              std::cout << " uvw " << geom_appro->getUVW() << " xyz " << geom_appro->getXYZ() << std::endl;
              std::cout << " e1 " << e1 << std::endl;
              std::cout << " e2 " << e2 << std::endl;
              std::cout << " nr  " << nr << std::endl;
	    }

      result = nr * grad_fr;
      if (debug) std::cout << " grad fr smooth case " << result << std::endl;
	  return;
	}
#endif
      // interpolated case
      result = fr_nodal_field.getGrad(e_appro);
      if (debug) std::cout << " grad fr interpolated case " << result << std::endl;
      return;
   }

  private:
   const lCrack& crack;
   const xfem::xLevelSet& lsn;
   const xfem::xLevelSet& lst;
   const double& rho;
   const xfem::xRegion& domain_for_integral;
   const xfem::xRegion& domain_for_integral_outer_layer;
   xfem::xLevelSet fr_nodal_field;
};

class xcEvalQGlobalVector : public xfem::xEval<xtensor::xVector<>>
{
  public:
   xcEvalQGlobalVector(const xfem::xVectorField& d, const xfem::xEval<double>& f) : eval_q_dir(d), eval_fr(f) {}

   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& result) const override
   {
      const bool debug = false;
      mEntity* e_appro = geo_appro->getEntity();
      xtensor::xPoint uvw_appro = geo_appro->getUVW();
      double fr;
      eval_fr(geo_appro, geo_integ, fr);
      result = eval_q_dir.getVal(e_appro, uvw_appro);
      if (debug) std::cout << " xyz for q_dir " << geo_integ->getXYZ() << std::endl;
      if (debug) std::cout << "q_dir " << result << std::endl;
      result *= fr;
      return;
   }

  private:
   const xfem::xVectorField& eval_q_dir;
   const xfem::xEval<double>& eval_fr;
};

class xcEvalGradQGlobalVector : public xfem::xEval<xtensor::xTensor2<>>
{
  public:
   xcEvalGradQGlobalVector(const xfem::xVectorField& d, const xfem::xEval<double>& f, const xfem::xEval<xtensor::xVector<>>& fg)
       : eval_q_dir(d), eval_fr(f), eval_grad_fr(fg)
   {
   }

   void operator()(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ, result_type& result) const override
   {
      const bool debug = false;
      mEntity* e_appro = geo_appro->getEntity();
      xtensor::xPoint uvw_appro = geo_appro->getUVW();
      // le terme ci-dessous est négligé dans Parks
      result = eval_q_dir.getGrad(e_appro);

      if (debug) std::cout << " xyz for grad_q_dir " << geo_integ->getXYZ() << std::endl;
      if (debug) std::cout << "grad_q_dir " << result << std::endl;

      double fr;
      eval_fr(geo_appro, geo_integ, fr);
      xtensor::xVector<> grad_fr;
      eval_grad_fr(geo_appro, geo_integ, grad_fr);

      result *= fr;

      result += tensor_product(eval_q_dir.getVal(e_appro, uvw_appro), grad_fr);
      return;
   }

  private:
   const xfem::xVectorField& eval_q_dir;
   const xfem::xEval<double>& eval_fr;
   const xfem::xEval<xtensor::xVector<>>& eval_grad_fr;
};

}  // namespace xcrack

#endif
