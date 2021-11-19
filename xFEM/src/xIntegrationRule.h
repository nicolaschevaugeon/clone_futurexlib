/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __INTEGRATORS_ABST_H
#define __INTEGRATORS_ABST_H
// std
#include <iostream>

// xmapping
#include "xReferenceElement.h"
// aomd
#include "AOMD.h"
#include "mEntity.h"

// xfem
#include "xGeomElem.h"
#include "xPartition.h"

// xinterface aomd
#include "xAttachedDataManagerAOMD.h"

namespace xfem
{
class xCommandOnGeomElem;
class xLevelSet;

class xRectangleIntegrator : public xquadrature::xIntegrator
{
  public:
   xRectangleIntegrator(AOMD::mEntity* bla) {}
   ~xRectangleIntegrator() override = default;
   int nbIntegrationPoints(int level) const override;
   void iPoint(int i, int level, double& u, double& v, double& w, double& weight) const override;
   static const int levelmax = 8;
};

// Two integrations : Gauss method and Rectangle method
class xIntegrationRule
{
  public:
   enum integrationRuleType
   {
      GAUSS,
      RECTANGLE
   };
   virtual ~xIntegrationRule()
   {
      if (integrator) delete integrator;
   };
   xIntegrationRule();
   virtual void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const = 0;
   void setIntegrator(integrationRuleType rule) { currentRule = rule; };
   xquadrature::xIntegrator* getIntegrator(AOMD::mEntity* e_integ) const
   {
      if (integrator)
      {
         delete integrator;
         integrator = nullptr;
      }
      switch (currentRule)
      {
         case GAUSS:
         {
            integrator = new xquadrature::xGaussIntegrator(e_integ);
            break;
         }
         case RECTANGLE:
         {
            integrator = new xRectangleIntegrator(e_integ);
            break;
         }
      }
      return integrator;
   };

  protected:
   integrationRuleType currentRule;
   mutable xquadrature::xIntegrator* integrator;
};

// gives local coordinates of integration points, for integration with rectangle method
std::vector<xtensor::xPoint> integpoint(int level);
void integpointnewlevelG(int u0, int v0, int* u, int* v);
void integpointnewlevelD(int u0, int v0, int* u, int* v);

//////////

class xIntegrationRuleBasic : public xIntegrationRule
{
  public:
   xIntegrationRuleBasic(int d = 0) : degree(d), nbpoint(0), filter(xAcceptAll()) {}
   template <class T>
   xIntegrationRuleBasic(const T& f, int d = 0) : degree(d), nbpoint(0), filter(f)
   {
   }
   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;
   int getNbPoints() const { return nbpoint; }
   void resetPointCounter() { nbpoint = 0; }

  private:
   int degree;
   mutable int nbpoint;
   xEntityFilter filter;
};

class xIntegrationRuleNodal : public xIntegrationRule
{
  public:
   xIntegrationRuleNodal() : filter(xAcceptAll()) {}
   template <class T>
   xIntegrationRuleNodal(const T& f) : filter(f)
   {
   }
   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

  private:
   xEntityFilter filter;
};

class xIntegrationRulePartition : public xIntegrationRule
{
  public:
   xIntegrationRulePartition(xGetPartition _getpartition, int de = 0)
       : degree(de), nbpoint(0), filter(xAcceptAll()), getpartition(_getpartition)
   {
   }
   xIntegrationRulePartition(xGetPartition _getpartition, const xEntityFilter& f, int de = 0)
       : degree(de), nbpoint(0), filter(f), getpartition(_getpartition)
   {
   }
   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;
   /// return the total number of gauss points used at the last integration
   /*! done by calling accept.
     (Return the sum of the number of gauss point used in each partition) */
   int getNbPoints() const { return nbpoint; }
   void resetPointCounter() { nbpoint = 0; }

  private:
   int degree;
   mutable int nbpoint;
   xEntityFilter filter;
   xGetPartition getpartition;
};

class xIntegrationRulePartitionBoundary : public xIntegrationRule
{
  public:
   xIntegrationRulePartitionBoundary(xGetPartition _getpartition, int de = 0)
       : getpartition(_getpartition), degree(de), filter(xAcceptAll())
   {
   }

   template <class T>
   xIntegrationRulePartitionBoundary(xGetPartition _getpartition, const T& f, int de = 0)
       : getpartition(_getpartition), degree(de), filter(f)
   {
   }

   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

  private:
   xGetPartition getpartition;
   int degree;
   xEntityFilter filter;
};

class xIntegrationRulePartitionExport : public xIntegrationRule
{
  public:
   xIntegrationRulePartitionExport(xGetPartition _getpartition, const xLevelSet* lst_, int de = 0)
       : degree(de), filter(xAcceptAll()), lst(lst_), getpartition(_getpartition)
   {
   }

   template <class T>
   xIntegrationRulePartitionExport(xGetPartition _getpartition, const T& f, int de = 0)
       : degree(de), filter(f), getpartition(_getpartition)
   {
   }

   void accept(xCommandOnGeomElem& command, AOMD::mEntity* e_integ) const override;

  private:
   int degree;
   xEntityFilter filter;
   xLevelSet const* lst;
   xGetPartition getpartition;
};

class xNodalIntegrationPoints : public xquadrature::xIntegrator
{
   // pEntity ent;
   AOMD::mEntity::mType entityType;

  public:
   xNodalIntegrationPoints(AOMD::mEntity*);
   int nbIntegrationPoints(int order) const override;
   void iPoint(int i, int order, double& u, double& v, double& w, double& weight) const override;
};

}  // namespace xfem

#endif
