/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef X_EXTEND_SHAPE_FCT_IMP_H
#define X_EXTEND_SHAPE_FCT_IMP_H
// xmapping
#include "xMappingBuilderHolder.h"
namespace xfem
{
template <typename DIST, typename VT>
template <typename ITERFRONT>
xExtendShapeGenerator<DIST, VT>::xExtendShapeGenerator(
    ITERFRONT beg, ITERFRONT end, xEntityToEntity upper_, xfem::xField<VT> &vgama_,
    std::function<void(AOMD::mEntity *, const xtensor::xPoint &, xmapping::xMapping *, xIntegrationRuleStoredDataManager &)>
        integration_order_function_object_,
    xfem::xIntegrationRuleStoredDataManager &gp_storage_)
    : distance_kd_tree(beg, end),
      upper(upper_),
      vgama(vgama_),
      gp_storage(gp_storage_),
      integration_order_function_object(integration_order_function_object_)
{
}

template <typename DIST, typename VT>
xExtendShapeFcts *xExtendShapeGenerator<DIST, VT>::generateExtendShapeFcts(AOMD::mEntity *e)
{
#ifndef NDEBUG
   // local
   // const bool debug_mode = xdebug_flag;
   const bool debug_mode = false;
#endif

   // generate xExtendShapeFcts and attach it to e
   xExtendShapeFcts &extention = data.setData(*e);

   // set geom elem
   xGeomElem elem(e);

   xmapping::xMapping *mapping = const_cast<xmapping::xMapping *>(elem.getMapping());

   // set gauss point and attach them to this element based on CDG
   integration_order_function_object(e, elem.getCDGxyz(), mapping, gp_storage);

   // set integrator based on integration_order_function_object generation
   xGaussPoints *agp = gp_storage.getStored(*e);
   if (!agp)
   {
      throw -78545;
   }
   xfem::xStoredIntegrationPoints integration_points(*agp);
   elem.setIntegrator(&integration_points);

   // get number of gauss point for e
   int dumy = 1;
   elem.SetIntegrationPointNumberForDegree(dumy);
   const int nb_gp = elem.GetNbIntegrationPoints();
#ifndef NDEBUG
   if (debug_mode)
   {
      std::cout << "Element " << e->getId() << " with CDG " << elem.getCDGxyz() << " got " << nb_gp << " gaus point\n";
   }
#endif

   // intermediate storage
   typedef std::multimap<AOMD::mEntity *, std::pair<int, xtensor::xPoint>> asso_t;
   typedef asso_t::iterator asso_it_t;
   asso_t asso_elem_gp;
   std::pair<asso_it_t, asso_it_t> equal_asso_elem_gp;
   std::set<AOMD::mEntity *> elem_set;

   // loop on gauss point to compute nearest stuff
   for (int i = 0; i < nb_gp; ++i)
   {
      // get gauss point coordinate
      elem.setUVW(i);
      auto gp_xyztmp = elem.getXYZ();
      xtensor::xPoint gp_xyz(gp_xyztmp(0), gp_xyztmp(1), gp_xyztmp(2));

      // get closest on front
      double distance;
      xtensor::xPoint nearest_point;
      AOMD::mEntity *nearest_entity;
      distance_kd_tree.nearestEntityDistance(gp_xyz, nearest_entity, nearest_point, distance);

      // element to work on
      AOMD::mEntity *nearest_true_entity = upper(nearest_entity);

      // store into set
      elem_set.insert(nearest_true_entity);

      // store into multimap
      asso_elem_gp.insert(make_pair(nearest_true_entity, make_pair(i, nearest_point)));

#ifndef NDEBUG
      if (debug_mode)
      {
         std::cout << "asso : key " << nearest_true_entity << " pair : gp id " << i << " nearest_point " << nearest_point << "\n";
      }
#endif
   }

   // fem description
   xfem::xFiniteElement FEM;

   // extention search variable
   std::pair<xExtendShapeFcts::extShapeFct_it_t, bool> gaussvalue;
   xExtendShapeFcts::extGaussVal_t dummy;
   std::pair<xExtendShapeFcts::extGaussVal_it_t, bool> inserted;
   double val;

   // loop on concerned entity to fill extention structure
   std::set<AOMD::mEntity *>::const_iterator it_elem = elem_set.begin();
   std::set<AOMD::mEntity *>::const_iterator it_elem_end = elem_set.end();
   for (; it_elem != it_elem_end; ++it_elem)
   {
      // nearest element
      AOMD::mEntity *ne = *it_elem;

      // get keys and function for this element
      FEM.setKeysAndFcts(ne, vgama.begin(), vgama.end());

      // set geom elem for this element
      xGeomElem nelem(ne);

      // get "gaus point" on that element
      equal_asso_elem_gp = asso_elem_gp.equal_range(ne);

      // loop on keys and fcts
      xfem::xFiniteElement::iterKey itk = FEM.beginKey();
      xfem::xFiniteElement::iterFct itf = FEM.beginFct();
      xfem::xFiniteElement::iterKey itk_end = FEM.endKey();
      for (; itk != itk_end; ++itk, ++itf)
      {
         // look if key allready present in the extention map and if not insert a empty
         gaussvalue = extention.associated_keys_and_value.insert(std::make_pair(*itk, dummy));

         // loop on "gaus point" on that element
         for (asso_it_t itgp = equal_asso_elem_gp.first; itgp != equal_asso_elem_gp.second; ++itgp)
         {
            // set "gauss point" position
            nelem.setUVWForXYZ(itgp->second.second);

            // compute value for this shape function for this "gauss point"
            (*itf)->getVal(&nelem, &nelem, val);

            // store in extention
            inserted = gaussvalue.first->second.insert(std::make_pair(itgp->second.first, val));

#ifndef NDEBUG
            if (debug_mode)
            {
               std::cout << "extend : new key " << *itk << " new pair : gp id " << itgp->second.first << " shape function value "
                         << val << " container size" << gaussvalue.first->second.size() << "\n";
            }
#endif
            // check
            if (!inserted.second)
            {
               throw -999;
            }
         }
      }
   }

   return &extention;
}

template <typename DIST, typename VT>
template <typename ITER>
void xExtendShapeGenerator<DIST, VT>::cleanExtendShapeFcts(ITER begin, ITER end)
{
   // loop on element
   for (; begin != end; ++begin)
   {
      AOMD::mEntity *e = *begin;

      // first treat element where ever it have a attached partition
      cleanExtendShapeFctsElem(e);

      // look for attached partition
      if (xfem::xMesh *m = xfem::xMesh::get_partition().getData(*e))
      {
         this->cleanExtendShapeFcts(m->begin(m->dim()), m->end(m->dim()));
      }
   }
   return;
}

template <typename DIST, typename VT>
void xExtendShapeGenerator<DIST, VT>::cleanExtendShapeFctsElem(AOMD::mEntity *e)
{
   gp_storage.delStored(*e);
   data.deleteData(*e);
   return;
}

}  // namespace xfem

#endif
