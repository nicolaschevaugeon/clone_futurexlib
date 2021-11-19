/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xApproxFunctionCorrectedShifted.h"

#include "xDebug.h"

using std::cerr;
using std::cout;
using std::endl;

namespace xfem
{
void xApproxFunctionEnrichedXFEMShiftedCorr::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double& res) const
{
   const bool debug = xdebug_flag;
   double b, e;
   base->getVal(geo_appro, geo_integ, b);
   //  cout<<geo_integ->getXYZ()<<endl;
   enrichment->getVal(geo_appro, geo_integ, e);

   xGeomElem geo_entity(enti);
   geo_entity.setUVW(uvwInterpo);
   xGeomElem geo_appro2(geo_appro->getEntity());
   geo_appro2.setUVWForXYZ(geo_entity.getXYZ());
   double ee;
   enrichment->getVal(&geo_appro2, &geo_entity, ee);
   if (enti) enti->print();

   double ramp = 0;
   (*ev_ramp)(geo_appro, geo_integ, ramp);

   res = b * (e - ee) * ramp;

   //  cout<<"xApproxFunctionEnrichedXFEMShiftedCorr\n";

   //  cout<<"GetVal base val "<<b<<" enric val "<<e<<" enrich shift "<<ee<<" shifted "<<e-ee<<" final "<<res<<endl;

   if (debug)
   {
      cout << " enrichment at point " << geo_appro->getXYZ() << " is " << e;
      cout << " ScalarxApproxFunctionEnrichedXFEM::getVal " << res << endl;
   }
}

void xApproxFunctionEnrichedXFEMShiftedCorr::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                     xtensor::xVector<>& res) const
{
   cout << "xApproxFunctionEnrichedXFEMShiftedCorr-------\n";
   const bool debug = xdebug_flag;
   xtensor::xVector<> bg, eg;
   double es, bs;

   base->getGrad(geo_appro, geo_integ, bg);
   enrichment->getVal(geo_appro, geo_integ, es);
   enrichment->getGrad(geo_appro, geo_integ, eg);
   base->getVal(geo_appro, geo_integ, bs);
   //  res  = bg * es + eg * bs;

   xGeomElem geo_entity(enti);
   geo_entity.setUVW(uvwInterpo);
   xGeomElem geo_appro2(geo_appro->getEntity());
   geo_appro2.setUVWForXYZ(geo_entity.getXYZ());
   double ee;
   enrichment->getVal(&geo_appro2, &geo_entity, ee);

   double ramp = 0;
   (*ev_ramp)(geo_appro, geo_integ, ramp);
   xtensor::xVector<> gramp;
   (*ev_gramp)(geo_appro, geo_integ, gramp);

   //  res  = bg * (es-ee) + eg * bs;
   res = bg * (es - ee) * ramp + eg * bs * ramp + gramp * (es - ee) * bs;

   //  cout<<"GetVal base val "<<bs<<" enric val "<<es<<" enrich shift "<<ee<<" shifted "<<e-ee<<" final "<<res<<endl;

   // doit etre corrige pour fct en shifted
   //  throw;

   if (debug)
   {
      cout << " enrichment  at point " << geo_appro->getXYZ() << " is " << es << endl;
      cout << " grad enrichment " << eg << endl;
      cout << " ScalarxApproxFunctionEnrichedXFEM::getGrad " << res << endl;
   }

   return;
}

void xApproxFunctionEnrichedXFEMShiftedCorr::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                    xtensor::xVector<>& res) const
{
   xtensor::xVector<> bg;
   double es;
   base->getVal(geo_appro, geo_integ, bg);
   enrichment->getVal(geo_appro, geo_integ, es);

   xGeomElem geo_entity(enti);
   geo_entity.setUVW(uvwInterpo);
   xGeomElem geo_appro2(geo_appro->getEntity());
   geo_appro2.setUVWForXYZ(geo_entity.getXYZ());
   double ee;
   //  enti->print();
   //  cout<<geo_entity.getUVW()<<endl;
   //  if(enti->getLevel()==1) throw;

   //  cout<<"uvwInterpo "<<uvwInterpo(0)<<" "<<uvwInterpo(1)<<endl;
   enrichment->getVal(&geo_appro2, &geo_entity, ee);

   double ramp = 0;
   (*ev_ramp)(geo_appro, geo_integ, ramp);

   res = bg * (es - ee) * ramp;
   return;
}

void xApproxFunctionEnrichedXFEMShiftedCorr::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                     xtensor::xTensor2<>& res) const
{
   const bool debug = xdebug_flag;
   xtensor::xTensor2<> bt;
   xtensor::xVector<> eg, bv;
   double es;
   base->getGrad(geo_appro, geo_integ, bt);
   enrichment->getVal(geo_appro, geo_integ, es);
   enrichment->getGrad(geo_appro, geo_integ, eg);
   base->getVal(geo_appro, geo_integ, bv);

   xGeomElem geo_entity(enti);
   geo_entity.setUVW(uvwInterpo);
   xGeomElem geo_appro2(geo_appro->getEntity());
   geo_appro2.setUVWForXYZ(geo_entity.getXYZ());
   double ee;
   enrichment->getVal(&geo_appro2, &geo_entity, ee);

   double ramp = 0;
   (*ev_ramp)(geo_appro, geo_integ, ramp);
   xtensor::xVector<> gramp;
   (*ev_gramp)(geo_appro, geo_integ, gramp);

   //  res =bt  * (es-ee) + tensor_product(bv,eg);
   res = bt * (es - ee) * ramp + tensor_product(bv, eg * ramp) + tensor_product(bv, gramp * (es - ee));

   if (debug)
   {
      cout << " enrichment  at point " << geo_appro->getXYZ() << " is " << es << endl;
      cout << " grad enrichment " << eg << endl;
      cout << "val classical " << endl;
      cout << bv << endl;
      cout << "grad classical " << endl;
      cout << bt << endl;

      cout << "first part " << endl;
      cout << bt * es << endl;
      cout << "second part " << endl;
      cout << xtensor::tensor_product(eg, bv) << endl;
   }
}

}  // namespace xfem
