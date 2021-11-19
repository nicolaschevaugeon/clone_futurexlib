/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/
#ifndef XAPPROXFUNCTIONDERIVDISCOCTREE_
#define XAPPROXFUNCTIONDERIVDISCOCTREE_

#include <cassert>
#include <iostream>
// xmapping
#include "xLagrangeMapping.h"
// trellis
#include "mEdge.h"
#include "mVertex.h"

// xtensor
#include "xVector.h"

// xinterface
#include "xAttachedDataManagerAOMD.h"
// xfem
#include "xApproxFctPtr.h"
#include "xApproxFunction.h"
#include "xElement.h"
#include "xEntityFilter.h"
#include "xEvalStorage.h"
#include "xGeomElem.h"
#include "xMesh.h"
#include "xValKey.h"

namespace xfem
{
// Enrichment function for material interface if octree refinement
template <template <class> class DATAMANAGER>
class xScalarFunctionDerivDiscXFEMOctree : public xApproxFunction, private xEvalStorage
{
  public:
   xScalarFunctionDerivDiscXFEMOctree(const DATAMANAGER<AOMD::mEntity*>& was_created_by_,
                                      const DATAMANAGER<double>& levelset_value_, const DATAMANAGER<int>& cut_edge_,
                                      const DATAMANAGER<xMesh>& refined_elements_,
                                      const DATAMANAGER<AOMD::mEntity*>& root_entity_);
   void getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ, double&) const override;
   void getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>&) const override;

  private:
   void computeEnr(AOMD::mEntity* e, AOMD::mEntity* e_int, AOMD::mEntity* e_meshr, const xtensor::xPoint& vi, double& enr) const;
   const DATAMANAGER<AOMD::mEntity*>& was_created_by;
   const DATAMANAGER<double>& levelset_value;
   const DATAMANAGER<int>& cut_edge;
   const DATAMANAGER<xMesh>& refined_elements;
   const DATAMANAGER<AOMD::mEntity*>& root_entity;
   mutable xinterface::aomd::xAttachedDataManagerAOMD<std::vector<double>> eint2enrichment;
};

template <template <class> class DATAMANAGER>
xScalarFunctionDerivDiscXFEMOctree<DATAMANAGER>::xScalarFunctionDerivDiscXFEMOctree(
    const DATAMANAGER<AOMD::mEntity*>& was_created_by_, const DATAMANAGER<double>& levelset_value_,
    const DATAMANAGER<int>& cut_edge_, const DATAMANAGER<xMesh>& refined_elements_,
    const DATAMANAGER<AOMD::mEntity*>& root_entity_)
    : was_created_by(was_created_by_),
      levelset_value(levelset_value_),
      cut_edge(cut_edge_),
      refined_elements(refined_elements_),
      root_entity(root_entity_)
{
}

template <template <class> class DATAMANAGER>
void xScalarFunctionDerivDiscXFEMOctree<DATAMANAGER>::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                             double& res) const
{
   AOMD::mEntity* e = geo_appro->getEntity();
   AOMD::mEntity* e_int = geo_integ->getEntity();
   AOMD::mEntity* e_meshr = was_created_by.at(*e_int);

   const std::vector<double>* penri = eint2enrichment.getData(*e_int);
   if (!penri)
   {
      const size_t nnode = size_t(e_int->size(0));
      std::vector<double> tmp(nnode);  // computation of the enrichment function at 3 nodes of e_int
      for (size_t i = 0; i < nnode; i++)
      {
         xtensor::xPoint vi = static_cast<AOMD::mVertex*>(e_int->get(0, int(i)))->point();
         computeEnr(e, e_int, e_meshr, vi, tmp[i]);
      }
      eint2enrichment.setData(*e_int) = tmp;
      penri = eint2enrichment.getData(*e_int);
   }
   // cout << "  enri is  :  " <<  (*enri)[0] <<" ,  "<<  (*enri)[1] <<" ,  "<<  (*enri)[2] << endl;
   xfem::xElement elem_int(e_int);
   elem_int.setUvw(geo_integ->getUVW());
   res = elem_int.getInterpoSca(*penri);  // linear interpolation at Gauss point of elem_int
                                          // cout << " res is  : " << res << endl;
}

template <template <class> class DATAMANAGER>
void xScalarFunctionDerivDiscXFEMOctree<DATAMANAGER>::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                              xtensor::xVector<>& res) const
{
   // const bool debug = xdebug_flag;
   if (exist(geo_appro, geo_integ, res)) return;
   AOMD::mEntity* e = geo_appro->getEntity();
   AOMD::mEntity* e_int = geo_integ->getEntity();
   AOMD::mEntity* e_meshr = was_created_by.at(*e_int);

   const std::vector<double>* penri = eint2enrichment.getData(*e_int);
   if (!penri)
   {
      const size_t nnode = size_t(e_int->size(0));
      std::vector<double> tmp(nnode);  // computation of the enrichment function at 3 nodes of e_int
      for (size_t i = 0; i < nnode; i++)
      {
         xtensor::xPoint vi = static_cast<AOMD::mVertex*>(e_int->get(0, int(i)))->point();
         computeEnr(e, e_int, e_meshr, vi, tmp[i]);
      }
      eint2enrichment.setData(*e_int) = tmp;
      penri = eint2enrichment.getData(*e_int);
   }
   xElement elem_int(e_int);
   elem_int.setUvw(geo_integ->getUVW());

   res = elem_int.getGradInterpoSca(*penri);
   return;
}

// Enrichment functions when octree refinement of the geometrical grid
template <template <class> class DATAMANAGER>
void xScalarFunctionDerivDiscXFEMOctree<DATAMANAGER>::computeEnr(AOMD::mEntity* e, AOMD::mEntity* e_int, AOMD::mEntity* e_meshr,
                                                                 const xtensor::xPoint& vi, double& enr) const
{
   int nnode = e->size(0);

   xElement elem(e);
   elem.xyz2uvw(vi);
   std::vector<double> lambda;
   elem.getFF(lambda);

   xElement elemr(e_meshr);
   elemr.xyz2uvw(vi);
   std::vector<double> lambdar;
   elemr.getFF(lambdar);

   std::vector<double> vals(size_t(nnode), 0);
   for (int i = 0; i != nnode; ++i) vals[size_t(i)] = levelset_value.at(*static_cast<AOMD::mVertex*>(e->get(0, i)));
   std::vector<double> valsr(size_t(e_meshr->size(0)), 0);
   for (int i = 0; i != e_meshr->size(0); ++i)
      valsr[size_t(i)] = levelset_value.at(*static_cast<AOMD::mVertex*>(e_meshr->get(0, i)));

   enr = 0.;
   for (int i = 0; i < nnode; ++i) enr += fabs(vals[i]) * lambda[i];  //  Contribution of the mesh level set
   enr -= fabs(elemr.getInterpoSca(valsr));                           //  Contribution of the octree level set

   std::vector<AOMD::mEdge*> edge_not_cut(3, nullptr);  // if dim = 2

   std::vector<AOMD::mVertex*> vv(nnode, nullptr);
   for (int i = 0; i < nnode; ++i)
   {
      vv[i] = static_cast<AOMD::mVertex*>(e->get(0, i));
   }

   std::vector<AOMD::mEdge*> cut_edges(3, nullptr);
   for (int i = 0; i < 3; ++i)
   {
      AOMD::mEdge* edge = (AOMD::mEdge*)e->get(1, i);
      int ce = (cut_edge.getData(*edge)) ? cut_edge.at(*edge) : 0;
      if (ce == 1 || ce == 3)
      {
         //  1 if edge not cut, 2 if cut edge, 3 if ls crosses the edge only by 2 nodes
         edge_not_cut[i] = edge;
      }
   }

   for (int k = 0; k < e->size(1); k++)
   {
      if (edge_not_cut[k] != nullptr)
      {
         AOMD::mEdge* qnc = edge_not_cut[k];
         int i1 = k;
         int i2 = (k + 1) % nnode;

         enr += -(fabs(vals[i1]) * lambda[i1]) -
                (fabs(vals[i2]) *
                 lambda[i2]);  // we subtract the linear part to edge bubble function (corresponding to an edge not cut)

         double edgefunfac = lambda[i1] + lambda[i2];
         double chi = 0.;
         if (fabs(lambda[i1] + lambda[i2]) > 1.e-6) chi = lambda[i2] / (lambda[i1] + lambda[i2]);
         if (vv[i1] != qnc->get(0, 0)) chi = 1. - chi;
         chi = chi * 2 - 1.;  // coordinates edge : [-1;1]

         const xMesh* pMesh_refined = refined_elements.getData(*e);
         if (pMesh_refined)
         {
            const xMesh& Mesh_refined = refined_elements.at(*e);
            // xmapping::xLagrangeMapping mapEdgeTop(qnc);
            xmapping::xLagrangeMapping mapEdgeTop(xinterface::aomd::getRefElementType(*qnc), xinterface::aomd::getPoints(*qnc));
            double x = 0., y = 0., z = 0.;
            mapEdgeTop.eval(chi, 0., 0., x, y, z);
            // find the subedge on boundary of e
            std::vector<AOMD::mEdge*> qnc_subedges;
            for (AOMD::mEntity* pse : Mesh_refined.range(1))
               if (root_entity.at(*pse) == qnc) qnc_subedges.push_back(static_cast<AOMD::mEdge*>(pse));
            AOMD::mEdge* qnc_sub_target = nullptr;
            double u_qnc_sub_target = 0.;
            for (AOMD::mEdge* qnc_sub : qnc_subedges)
            {
               // xmapping::xLagrangeMapping mapEdgeCurrent(qnc_sub);
               xmapping::xLagrangeMapping mapEdgeCurrent(xinterface::aomd::getRefElementType(*qnc_sub),
                                                         xinterface::aomd::getPoints(*qnc_sub));
               double u, v, w;
               mapEdgeCurrent.invert(x, y, z, u, v, w);
               if (mapEdgeCurrent.inReferenceElement(u, v, w))
               {
                  qnc_sub_target = qnc_sub;
                  u_qnc_sub_target = u;
                  break;
               }
            }
            if (!qnc_sub_target) throw;
            // know compute the enrichment contrib.
            xElement elemq(qnc_sub_target);
            elemq.setUvw(xtensor::xPoint(u_qnc_sub_target, 0., 0.));
            AOMD::mVertex* vq0 = static_cast<AOMD::mVertex*>(qnc_sub_target->get(0, 0));
            AOMD::mVertex* vq1 = static_cast<AOMD::mVertex*>(qnc_sub_target->get(0, 1));
            std::vector<double> valq = {levelset_value.at(*vq0), levelset_value.at(*vq1)};
            const double b_edgei = elemq.getInterpoSca(valq);
            const double ab_edgei = fabs(b_edgei);
            enr += ab_edgei * edgefunfac;
         }
         else
         {
            xElement elemq(qnc);
            elemq.setUvw(xtensor::xPoint(chi, 0, 0));
            AOMD::mVertex* vq0 = static_cast<AOMD::mVertex*>(qnc->get(0, 0));
            AOMD::mVertex* vq1 = static_cast<AOMD::mVertex*>(qnc->get(0, 1));
            std::vector<double> valq(2, 0);
            valq[0] = levelset_value.at(*vq0);
            valq[1] = levelset_value.at(*vq1);
            const double b_edgei = elemq.getInterpoSca(valq);
            const double ab_edgei = fabs(b_edgei);
            enr += ab_edgei * edgefunfac;
         }
      }
   }
   return;
}

}  // namespace xfem

#endif
