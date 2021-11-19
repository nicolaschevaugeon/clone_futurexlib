/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XINTEGRATION_RULE_STORED_H
#define _XINTEGRATION_RULE_STORED_H

// xquadrature
#include "xIntegrator.h"
// xfem include
#include "xIntegrationRule.h"

namespace xfem
{
using xGaussPoints = std::vector<std::pair<xtensor::xPoint, double>>;

/// Class to get the integration points attached to the elements
class xStoredIntegrationPoints : public xquadrature::xIntegrator
{
   const xGaussPoints& gauss_points;

  public:
   xStoredIntegrationPoints(const xGaussPoints& gauss_points) : gauss_points(gauss_points) {}
   int nbIntegrationPoints(int order) const override { return gauss_points.size(); }
   void iPoint(int i, int order, double& u, double& v, double& w, double& weight) const override
   {
      const xtensor::xPoint& uvw = gauss_points[i].first;
      u = uvw(0);
      v = uvw(1);
      w = uvw(2);
      weight = gauss_points[i].second;
   }
};

class xIntegrationRuleStoredDataManager
{
  public:
   xGaussPoints* getStored(AOMD::mEntity& e);
   void delStored(AOMD::mEntity& e);
   xGaussPoints& createEmptyStorageOrGet(AOMD::mEntity& e);
   void clearStored();

  protected:
   xinterface::aomd::xAttachedDataManagerAOMD<xGaussPoints> gauss_points;
};

class xIntegrationRuleStoredBasic : public xIntegrationRule, public xIntegrationRuleStoredDataManager
{
  public:
   template <class T>
   xIntegrationRuleStoredBasic(const T& f, int d = 0) : degree(d), filter(f), nbpoint(0)
   {
   }
   xIntegrationRuleStoredBasic(int d = 0) : xIntegrationRuleStoredBasic(xAcceptAll(), d) {}
   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

  private:
   const int degree;
   xEntityFilter filter;
   mutable int nbpoint;
};

/// Integration rule that uses stored integration points if they are present, and partition otherwise
class xIntegrationRuleStoredBasicFallBackPartition : public xIntegrationRule, public xIntegrationRuleStoredDataManager
{
  public:
   template <class T>
   xIntegrationRuleStoredBasicFallBackPartition(xGetPartition _getpartition, const T& f, int d = 0)
       : degree(d), filter(f), nbpoint(0), getpartition(_getpartition)
   {
   }
   xIntegrationRuleStoredBasicFallBackPartition(xGetPartition _getpartition, int d = 0)
       : xIntegrationRuleStoredBasicFallBackPartition(_getpartition, xAcceptAll(), d)
   {
   }
   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

  private:
   const int degree;
   xEntityFilter filter;
   mutable int nbpoint;
   xGetPartition getpartition;
};

class xIntegrationRuleStoredPartition : public xIntegrationRule, public xIntegrationRuleStoredDataManager
{
  public:
   template <class T>
   xIntegrationRuleStoredPartition(xGetPartition _getpartition, const T& f, int d = 0)
       : degree(d), filter(f), getpartition(_getpartition), nbpoint(0)
   {
   }
   xIntegrationRuleStoredPartition(xGetPartition _getpartition, int d = 0)
       : xIntegrationRuleStoredPartition(_getpartition, xAcceptAll(), d)
   {
   }
   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

  private:
   const int degree;
   xEntityFilter filter;
   xGetPartition getpartition;
   mutable int nbpoint;
};

///// Integration rule based on gauss points attached on the elements.
// class xIntegrationRuleStored : public xIntegrationRule {
//   public:
//     xIntegrationRuleStored(xLevelSet &ls_, int d = 0) : degree(d), filter(xAcceptAll()), ls(ls_) {}
//     template<class T>
//     xIntegrationRuleStored(xLevelSet &ls_, const T& f, int d = 0) : degree(d), filter(f), ls(ls_) {}
//     void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const
//     {
//       if (!filter(e_integ)) {return;}
//       Trellis_Util::LagrangeMappingBuilder builder;
//       *  mapping  = builder.BuildMapping(e_integ);
//       xStoredIntegrationPoints integration_points(e_integ);
//       xGeomElem geo_integ(e_integ);
//       if(!integration_points.isClassical()){
//         integration_points.filter(ls);
//         geo_integ.setMapping(mapping);
//         geo_integ.setIntegrator(&integration_points);
// //         xGeomElem geo_integ2(e_integ, mapping, &integration_points);
// //         geo_integ=geo_integ2;
//       }
//       geo_integ.SetIntegrationPointNumberForDegree(degree);///degree not relevant if stored...
//       command.execute(&geo_integ);
//       delete mapping;
//     }
//
//
//
//   private:
//     xEntityFilter filter;
//     int degree;
//     xLevelSet &ls;
// };

}  // namespace xfem

#endif
