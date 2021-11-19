/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef __COMMAND_ongeomelem_H
#define __COMMAND_ongeomelem_H

#include "mAttachableDataContainer.h"  ///
#include "xGeomElem.h"
#include "xIntegrator.h"
//#include "mAOMD.h" ///

namespace xfem
{
class xForm;

class xCommandOnGeomElem
{
  public:
   virtual ~xCommandOnGeomElem() = default;
   virtual void openApproxElem(xGeomElem* g_appro);
   virtual void setIntegElem(xGeomElem* g_integ);
   virtual void execute(xGeomElem* g_integ) = 0;
   virtual void closeApproxElem(xGeomElem* g_appro);
   xGeomElem* getApproxElem() { return geom_appro; };

  protected:
   xGeomElem* geom_appro;
   xGeomElem* geom_integ;
};

template <typename Eval, typename Compare = std::less<typename Eval::result_type>>
class xEvalMaxCommand : public xfem::xCommandOnGeomElem
{
  public:
   xEvalMaxCommand(Eval& e, typename Eval::result_type& r) : eval(e), res(r) {}
   void execute(xGeomElem* geo_integ) override
   {
      const int nb = geo_integ->GetNbIntegrationPoints();
      for (int k = 0; k < nb; k++)
      {
         geo_integ->setUVW(k);
         if (geom_appro->getEntity() != geo_integ->getEntity())
            geom_appro->setUVWForXYZ(geo_integ->getXYZ());
         else
            geom_appro->setUVW(geo_integ->getUVW());
         eval(geom_appro, geo_integ, val);
         res = max(val, res, op);
      }
   }

  private:
   Eval& eval;
   typename Eval::result_type val, &res;
   Compare op;
};

class xIntegrateFormCommand : public xCommandOnGeomElem
{
  public:
   xIntegrateFormCommand(xForm* f);
   void execute(xGeomElem* geo_integ) override;
   xForm* getForm() { return form; }

  protected:
   xForm* form;
};

template <class Eval>
class xIntegrateEvalCommand : public xCommandOnGeomElem
{
  public:
   xIntegrateEvalCommand(const Eval& e, typename Eval::result_type& r) : eval(e), res(r) {}

   void execute(xGeomElem* geo_integ) override
   {
      const int nb = geo_integ->GetNbIntegrationPoints();
      // cout << "EXECUTE " << nb << endl;
      for (int k = 0; k < nb; k++)
      {
         geo_integ->setUVW(k);
         if (geom_appro->getEntity() != geo_integ->getEntity())
            geom_appro->setUVWForXYZ(geo_integ->getXYZ());
         else
            geom_appro->setUVW(geo_integ->getUVW());
         eval(geom_appro, geo_integ, val);
         // cout << geo_integ->GetWeight() << " "<< geo_integ->GetDetJac() << endl;
         res += val * geo_integ->GetWeight() * geo_integ->GetDetJac();
      }
   }
   inline void reset_result() { res = typename Eval::result_type(); };
   inline const typename Eval::result_type& get_result() const { return res; };

  public:
   const Eval& getEvalFunctor() const { return eval; };

  private:
   const Eval& eval;
   typename Eval::result_type val, &res;
};

}  // namespace xfem

#endif
