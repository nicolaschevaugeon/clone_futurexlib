/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMUPDATER_H
#define _FMUPDATER_H

#include "FMUtil.h"

#ifdef HAVE_XGRAPH
// xgraph
#include "nodeAndConnectedEdge.h"
#endif

namespace xfastmarching
{
namespace internal
{
// DATAMANAGER represent the type of the data manager concept to use in this class
// MESHINTERFACE interface represent the type of the interface to mesh
// VECT  is the type used to represent vectors
// SCAL represent the level set arithmetic and consecuentely the FM aritmetic
// TRANSPORTED  represent the type od data transforted on the fly by the updater
template <template <class> class DATAMANAGER, class MESHINTERFACE, class VECT, class SCAL, class TRANSPORTED>
class FMupdaterGeneric
{
  public:
   typedef SCAL scal;
   typedef VECT vect;
   typedef TRANSPORTED transported;
   typedef MESHINTERFACE meshinterface_t;

   typedef typename MESHINTERFACE::vertex vertex;
   typedef typename MESHINTERFACE::edge edge;
   typedef typename MESHINTERFACE::face face;
   typedef typename MESHINTERFACE::region region;

   typedef std::vector<const region *> r_container_t;
   typedef std::vector<const face *> f_container_t;
   typedef DATAMANAGER<SCAL> entitystorage_t;
   typedef DATAMANAGER<VECT> gls_t;
   typedef DATAMANAGER<TRANSPORTED> vn_t;

   FMupdaterGeneric(const MESHINTERFACE &_mi, entitystorage_t &_ls, std::function<scal(const vertex &)> _Ffunc,
                    scal _epsilon_ratio)
       : status(_mi), mi(_mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(nullptr), pvn(nullptr)
   {
   }
   FMupdaterGeneric(const MESHINTERFACE &_mi, entitystorage_t &_ls, std::function<scal(const vertex &)> _Ffunc,
                    scal _epsilon_ratio, gls_t &_gls)
       : status(_mi), mi(_mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(&_gls), pvn(nullptr)
   {
   }
   FMupdaterGeneric(const MESHINTERFACE &_mi, entitystorage_t &_ls, std::function<scal(const vertex &)> _Ffunc,
                    scal _epsilon_ratio, gls_t &_gls, vn_t &_vn)
       : status(_mi), mi(_mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(&_gls), pvn(&_vn)
   {
   }

   void setLs(const vertex &v, const SCAL &p) { ls.setData(v) = p; }
   bool getLs(const vertex &v, SCAL &p) const
   {
      const SCAL *pp = ls.getData(v);
      if (pp)
      {
         p = *pp;
         return true;
      }
      else
         return false;
   }
   const SCAL &getLs(const vertex &v) const
   {
      const SCAL *pp = ls.getData(v);
      assert(pp);
      return *pp;
   }
   bool updatefg(const vertex &vx, const edge &e)
   {
      const SCAL Fx = Ffunc(vx);
      const vertex *vb = getothernodes(mi, vx, e);
      const VECT b = getcoord(*vb) - getcoord(vx);
      SCAL epsilon = xtool::xDataType<SCAL>::zero();
      VECT Grad;
      computeTrial(SCAL(xtool::xDataType<SCAL>::zero()), dual(b), Fx, epsilon, Grad);
      epsilon *= epsilon_ratio;
      return updatefp(vx, e, epsilon);
   }

   bool updatef(const vertex &vx, const edge &e) { return updatefp(vx, e, xtool::xDataType<SCAL>::zero()); }

   const scal L = std::numeric_limits<scal>::max();
   Status<MESHINTERFACE, vertex> status;

  protected:
   const MESHINTERFACE &mi;

   DATAMANAGER<SCAL> &ls;
   std::function<SCAL(const vertex &)> Ffunc;
   const SCAL epsilon_ratio;

   gls_t *pgls;
   vn_t *pvn;

   vect getcoord(const vertex &v) const
   {
      vect res;
      getCoord(mi, v, res);
      return res;
   }

   virtual bool updatefp(const vertex &vx, const edge &e, const scal &epsilon)
   {
      SCAL valx = getLs(vx);
      VECT px = getcoord(vx);
      SCAL Fx = Ffunc(vx);
      SCAL Trial;
      VECT GradTrial;
      VECT Grad;
      TRANSPORTED Trans;
      bool updated = false;
      const std::function<bool(const vertex *v)> is_known = [this](const vertex *v) { return (status(v) == FMStatus::K); };
      std::function<bool(double, double)> less = [epsilon](double v, double w) { return v < w + epsilon; };

      r_container_t e2r;
      f_container_t e2f;
      e2r.reserve(50);
      e2f.reserve(50);

      getRegions(mi, e, std::back_inserter(e2r));
      for (auto r : e2r)
      {
         const auto va = getothernodes(mi, vx, *r);
         if (std::all_of(va.begin(), va.end(), is_known))
         {
            const VECT T = {getLs(*va[0]), getLs(*va[1]), getLs(*va[2])};
            if (less(T(0), valx) || less(T(1), valx) || less(T(2), valx))
            {
               const std::array<VECT, 3> a = {getcoord(*va[0]) - px, getcoord(*va[1]) - px, getcoord(*va[2]) - px};
               const tensor3d<SCAL> dualG = dual(tensor3d<SCAL>{a[0], a[1], a[2]});
               if (computeTrial(T, dualG, Fx, Trial, GradTrial))
               {
                  if (less(Trial, valx))
                  {
                     updated = true;
                     FM_OUTPUT(" update region for " << vx.getId() << " with " << valx << " by " << Trial << " ("
                                                     << va[0]->getId() << "," << va[1]->getId() << "," << va[2]->getId() << ")"
                                                     << endl);
                     valx = Trial;
                     Grad = GradTrial;
                     if (pvn)
                     {
                        const std::array<TRANSPORTED, 3> trans_other{*(pvn->getData(*va[0])), *(pvn->getData(*va[1])),
                                                                     *(pvn->getData(*va[2]))};
                        auto gc = Grad * dualG;
                        Trans =
                            (trans_other[0] * gc(0) + trans_other[1] * gc(1) + trans_other[2] * gc(2)) / (gc(0) + gc(1) + gc(2));
                     }
                  }
               }
            }
         }
      }

      getFaces(mi, e, std::back_inserter(e2f));
      for (auto f : e2f)
      {
         const auto va = getothernodes(mi, vx, *f);
         if (std::all_of(va.begin(), va.end(), is_known))
         {
            const vector2d<SCAL> T = {getLs(*va[0]), getLs(*va[1])};
            if (less(T(0), valx) || less(T(1), valx))
            {
               const std::array<VECT, 2> a{getcoord(*va[0]) - px, getcoord(*va[1]) - px};
               const auto dualG = dual(tensor3d2d<SCAL>(a[0], a[1]));
               if (computeTrial(T, dualG, Fx, Trial, GradTrial))
               {
                  if (less(Trial, valx))
                  {
                     updated = true;
                     FM_OUTPUT(" update face for " << vx.getId() << " with " << valx << " by " << Trial << " (" << va[0]->getId()
                                                   << "," << va[1]->getId() << ")" << endl);
                     valx = Trial;
                     Grad = GradTrial;
                     if (pvn)
                     {
                        const std::array<TRANSPORTED, 2> trans_other{*(pvn->getData(*va[0])), *(pvn->getData(*va[1]))};
                        auto gc = Grad * dualG;
                        Trans = (trans_other[0] * gc(0) + trans_other[1] * gc(1)) / (gc(0) + gc(1));
                     }
                  }
               }
            }
         }
      }

      const vertex *vb = getothernodes(mi, vx, e);
      if (is_known(vb))
      {
         const SCAL T = getLs(*vb);
         if (less(T, valx))
         {
            const VECT b = getcoord(*vb) - px;
            if (computeTrial(T, dual(b), Fx, Trial, GradTrial))
            {
               if (less(Trial, valx))
               {
                  updated = true;
                  FM_OUTPUT(" update edge for " << vx.getId() << " with " << valx << " by " << Trial << " (" << vb->getId() << ")"
                                                << endl);
                  valx = Trial;
                  Grad = GradTrial;
                  if (pvn)
                  {
                     Trans = *(pvn->getData(*vb));
                  }
               }
            }
         }
      }

      if (updated)
      {
         setLs(vx, valx);
         if (pgls) pgls->setData(vx) = Grad;
         if (pvn) pvn->setData(vx) = Trans;
         return true;
      }
      return false;
   }
};

#ifdef HAVE_XGRAPH
// DATAMANAGER represent the type of the data manager concept to use in this class
// MESHINTERFACE interface represent the type of the interface to mesh
// VECT  is the type used to represent vectors
// SCAL represent the level set arithmetic and consecuentely the FM aritmetic
template <template <class> class DATAMANAGER, class MESHINTERFACE, class VECT, class SCAL>
class FMupdaterFlux
{
  public:
   typedef SCAL scal;
   typedef VECT vect;

   typedef typename MESHINTERFACE::vertex vertex;
   typedef typename MESHINTERFACE::edge edge;
   typedef typename MESHINTERFACE::face face;
   typedef typename MESHINTERFACE::region region;

   typedef std::vector<const region *> r_container_t;
   typedef std::vector<const face *> f_container_t;
   typedef DATAMANAGER<SCAL> entitystorage_t;
   typedef DATAMANAGER<VECT> gls_t;
   typedef xgraph::nodeFrom<vertex, SCAL, VECT, 3> node_from_t;
   typedef xgraph::nodeTo<vertex, SCAL, VECT> node_to_t;
   typedef DATAMANAGER<node_from_t> flux_from_t;
   typedef DATAMANAGER<node_to_t> flux_to_t;

   FMupdaterFlux(const MESHINTERFACE &_mi, entitystorage_t &_ls, std::function<scal(const vertex &)> _Ffunc, SCAL _epsilon_ratio,
                 flux_from_t &_flux_from, flux_to_t &_flux_to)
       : status(_mi), mi(_mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), flux_from(_flux_from), flux_to(_flux_to)
   {
   }

   void setLs(const vertex &v, const SCAL &p) { ls.setData(v) = p; }
   bool getLs(const vertex &v, SCAL &p) const
   {
      const SCAL *pp = ls.getData(v);
      if (pp)
      {
         p = *pp;
         return true;
      }
      else
         return false;
   }
   const SCAL &getLs(const vertex &v) const
   {
      const SCAL *pp = ls.getData(v);
      assert(pp);
      return *pp;
   }
   void setToKnown(const vertex *vk)
   {
      status(vk) = FMStatus::K;
      FM_OUTPUT(" set to known vertex " << vk->getId() << endl);
      node_from_t *from = flux_from.getData(*vk);
      if (from)
      {
         // loop on parent if any to set child in them
         auto gedges = from->getParentEdges();
         for (auto &gedge : gedges)
         {
            FM_OUTPUT(" update its parent " << gedge.first->getId() << endl);
            node_to_t &parent = flux_to.setData(*gedge.first);
            parent.insertChild(vk, gedge.second);
         }
         // store gradiant
         node_to_t &current = flux_to.setData(*vk);
         const VECT *grad = from->getAttachedData();
         if (grad)
         {
            current.setAttachedData(*grad);
            // to gain some memory
            // but info lost then
            from->delAttachedData();
         }
      }
   }
   bool updatefg(const vertex &vx, const edge &e)
   {
      const SCAL Fx = Ffunc(vx);
      const vertex *vb = getothernodes(mi, vx, e);
      const VECT b = getcoord(*vb) - getcoord(vx);
      SCAL epsilon = xtool::xDataType<SCAL>::zero();
      VECT Grad;
      computeTrial(SCAL(xtool::xDataType<SCAL>::zero()), dual(b), Fx, epsilon, Grad);
      epsilon *= epsilon_ratio;
      return updatefp(vx, e, epsilon);
   }

   bool updatef(const vertex &vx, const edge &e) { return updatefp(vx, e, xtool::xDataType<SCAL>::zero()); }

   const SCAL L = std::numeric_limits<SCAL>::max();
   Status<MESHINTERFACE, vertex> status;

  private:
   const MESHINTERFACE &mi;
   DATAMANAGER<SCAL> &ls;
   std::function<SCAL(const vertex &)> Ffunc;
   const SCAL epsilon_ratio;
   flux_from_t &flux_from;
   flux_to_t &flux_to;

   bool updatefp(const vertex &vx, const edge &e, const SCAL &epsilon)
   {
      const SCAL one = xtool::xDataType<SCAL>::one();
      SCAL valx = getLs(vx);
      VECT px = getcoord(vx);
      SCAL Fx = Ffunc(vx);
      SCAL Trial;
      VECT GradTrial;
      VECT Grad;
      node_from_t From;
      bool updated = false;
      auto is_known = [this](const vertex *v) {
         // return status(v);
         // return ( status(v) && status(v) != FMStatus::T );
         return (status(v) == FMStatus::K);
      };
      auto less = [epsilon](double v, double w) { return v < w + epsilon; };

      r_container_t e2r;
      f_container_t e2f;
      e2r.reserve(50);
      e2f.reserve(50);

      getRegions(mi, e, std::back_inserter(e2r));
      for (auto r : e2r)
      {
         const auto va = getothernodes(mi, vx, *r);
         if (std::all_of(va.begin(), va.end(), is_known))
         {
            const VECT T = {getLs(*va[0]), getLs(*va[1]), getLs(*va[2])};
            if (less(T(0), valx) || less(T(1), valx) || less(T(2), valx))
            {
               const std::array<VECT, 3> a = {getcoord(*va[0]) - px, getcoord(*va[1]) - px, getcoord(*va[2]) - px};
               const tensor3d<SCAL> dualG = dual(tensor3d<SCAL>{a[0], a[1], a[2]});
               if (computeTrial(T, dualG, Fx, Trial, GradTrial))
               {
                  if (less(Trial, valx))
                  {
                     updated = true;
                     FM_OUTPUT(" update region for " << vx.getId() << " with " << valx << " by " << Trial << " ("
                                                     << va[0]->getId() << "," << va[1]->getId() << "," << va[2]->getId() << ")"
                                                     << endl);
                     valx = Trial;
                     auto gc = GradTrial * dualG;
                     gc *= one / (gc(0) + gc(1) + gc(2));
                     From.set(va, &(gc(0)));
                     Grad = GradTrial;
                  }
               }
            }
         }
      }

      getFaces(mi, e, std::back_inserter(e2f));
      for (auto f : e2f)
      {
         const auto va = getothernodes(mi, vx, *f);
         if (std::all_of(va.begin(), va.end(), is_known))
         {
            const vector2d<SCAL> T = {getLs(*va[0]), getLs(*va[1])};
            if (less(T(0), valx) || less(T(1), valx))
            {
               const std::array<VECT, 2> a{getcoord(*va[0]) - px, getcoord(*va[1]) - px};
               const auto dualG = dual(tensor3d2d<SCAL>(a[0], a[1]));
               if (computeTrial(T, dualG, Fx, Trial, GradTrial))
               {
                  if (less(Trial, valx))
                  {
                     updated = true;
                     FM_OUTPUT(" update face for " << vx.getId() << " with " << valx << " by " << Trial << " (" << va[0]->getId()
                                                   << "," << va[1]->getId() << ")" << endl);
                     valx = Trial;
                     auto gc = GradTrial * dualG;
                     gc *= one / (gc(0) + gc(1));
                     From.set(va, &(gc(0)));
                     Grad = GradTrial;
                  }
               }
            }
         }
      }

      const vertex *vb = getothernodes(mi, vx, e);
      if (is_known(vb))
      {
         const SCAL T = getLs(*vb);
         if (less(T, valx))
         {
            const VECT b = getcoord(*vb) - px;
            if (computeTrial(T, dual(b), Fx, Trial, GradTrial))
            {
               if (less(Trial, valx))
               {
                  updated = true;
                  FM_OUTPUT(" update edge for " << vx.getId() << " with " << valx << " by " << Trial << " (" << vb->getId() << ")"
                                                << endl);
                  valx = Trial;
                  From.set(vb, one);
                  Grad = GradTrial;
               }
            }
         }
      }

      if (updated)
      {
         setLs(vx, valx);
         auto &gfrom = flux_from.setData(vx);
         gfrom = From;
         gfrom.setAttachedData(Grad);
         return true;
      }
      return false;
   }

   VECT getcoord(const vertex &v) const
   {
      VECT res;
      getCoord(mi, v, res);
      return res;
   };
};
#endif

}  // namespace internal
}  // namespace xfastmarching
#endif
