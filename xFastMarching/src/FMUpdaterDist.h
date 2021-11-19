/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMUPDATERDIST_H
#define _FMUPDATERDIST_H

#include "FMUtil.h"
#include "linearalgebra3dPackUnPack.h"
#include "xMPIDataType.h"

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
class FMupdaterGenericDist
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
   typedef DATAMANAGER<TRANSPORTED> vn_t;
   typedef std::tuple<SCAL, VECT, TRANSPORTED> tuple_data_t;
   typedef std::unordered_map<int, tuple_data_t> remote_tuple_data_t;
   typedef std::tuple<bool, SCAL, remote_tuple_data_t> tuple_data_bnd_t;
   typedef DATAMANAGER<tuple_data_bnd_t> data_bnd_t;

   FMupdaterGenericDist(const MESHINTERFACE &_mi, entitystorage_t &_ls, std::function<scal(const vertex &)> _Ffunc,
                        SCAL _epsilon_ratio)
       : status(_mi), mi(_mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(nullptr), pvn(nullptr)
   {
   }
   FMupdaterGenericDist(const MESHINTERFACE &_mi, entitystorage_t &_ls, std::function<scal(const vertex &)> _Ffunc,
                        SCAL _epsilon_ratio, gls_t &_gls)
       : status(_mi), mi(_mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(&_gls), pvn(nullptr)
   {
   }
   FMupdaterGenericDist(const MESHINTERFACE &_mi, entitystorage_t &_ls, std::function<scal(const vertex &)> _Ffunc,
                        SCAL _epsilon_ratio, gls_t &_gls, vn_t &_vn)
       : status(_mi), mi(_mi), ls(_ls), Ffunc(_Ffunc), epsilon_ratio(_epsilon_ratio), pgls(&_gls), pvn(&_vn)
   {
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

   void pack(const vertex &v, xtool::xMpiInputBuffer &buff)
   {
      tuple_data_bnd_t *pdbound = data_bnd.getData(v);
      if (pdbound)
      {
         assert(std::get<0>(*pdbound));
         SCAL valx;
         getLs(v, valx);
         buff.pack(&valx, 1, xtool::xMPIDataType<scal>());
         if (pgls)
         {
            VECT *pgrad = pgls->getData(v);
            if (!pgrad)
            {
               if (status(&v) != FMStatus::KF) throw -34567;
            }
            ::pack(*pgrad, buff);
         }

         if (pvn)
         {
            TRANSPORTED *Trans = pvn->getData(v);
            if (!Trans)
            {
               if (status(&v) != FMStatus::KF) throw -34568;
            }
            ::pack(*Trans, buff);
         }
      }
      else
         throw -12679;
   }
   std::pair<scal, const vertex *> unPack(const xtool::xMpiOutputBuffer &buff, int receivedfrom)
   {
      const vertex *v;
      buff.unPack(&v, 1, MPI_AINT);
      FM_OUTPUT("unpack " << v->getId() << " from " << receivedfrom << endl);
      tuple_data_bnd_t *pdbound = data_bnd.getData(*v);
      if (pdbound)
      {
         tuple_data_bnd_t &dbound = *pdbound;
         remote_tuple_data_t &rd = std::get<2>(dbound);
         scal oldv = L;
         auto itf = rd.find(receivedfrom);
         if (itf != rd.end()) oldv = std::get<0>(itf->second);
         auto &d = rd[receivedfrom];
         buff.unPack(&std::get<0>(d), 1, xtool::xMPIDataType<SCAL>());
         if (pgls) ::unPack(std::get<1>(d), buff);
         if (pvn) ::unPack(std::get<2>(d), buff);
         // std::get < 2 >(d).unPack(buff);

         // if the received val is the same as the old one, nothing change and the minimum is the same
         FM_OUTPUT("old " << oldv << " new " << std::get<0>(d) << " last min " << std::get<1>(dbound) << endl);
         if (std::get<0>(d) == oldv) return std::make_pair(L, nullptr);

         // compute the min of current received value from remote proc.
         // ==========================================================
         // One may say that this operation is done to many time and should have been done when
         // all message are received and the map is full, but:
         // * Doing this afterward means having a way to retrieve boundary vertex that did receive something
         // and this may have a huge cost.
         // * Doing this here does not mean that operation will be done so many time. The probability for a shared node
         // to receive a front from many process at the same time is not so hight. First because nodes are for most of them
         // share in between 2 proc and then only receive one message. Second the node shared must be on the skeleton to receive
         // many front, otherwise only one front, likely from one proc will be received.
         scal min_valr = L;
         for (auto &id : rd)
         {
            if (std::get<0>(id.second) < min_valr) min_valr = std::get<0>(id.second);
         }
         std::get<1>(dbound) = min_valr;
         FM_OUTPUT(" new min " << std::get<1>(dbound) << endl);

         // only if received val is the min of all received so far this function
         // return a valid pair. Otherwise it give a pair with a null pointer to say
         // that this point is not to be added to the stack
         if (std::get<0>(d) > min_valr)
            return std::make_pair(L, nullptr);
         else
         {
            FM_OUTPUT(" select as new " << std::get<0>(d) << " <= " << min_valr << endl);

            std::get<0>(dbound) = false;
            if (pgls) pgls->setData(*v) = std::get<1>(d);
            if (pvn) pvn->setData(*v) = std::get<2>(d);
            return std::make_pair(min_valr, v);
         }
      }
      else
         throw -12679;
   }

   bool getGlobalLs(const vertex &v, SCAL &p) const
   {
      // security
      p = L;
      // get if exist
      bool g = getLs(v, p);
      // boundary
      const tuple_data_bnd_t *pdbound = data_bnd.getData(v);
      if (pdbound)
      {
         p = std::min(std::get<1>(*pdbound), p);
         g = true;
      }
      return g;
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
   void setLsKnown(const vertex &v, SCAL &p)
   {
      getLs(v, p);
      tuple_data_bnd_t *pdbound = data_bnd.getData(v);
      if (pdbound)
      {
         tuple_data_bnd_t &dbound = *pdbound;
         std::get<0>(dbound) = false;
         std::get<1>(dbound) = p;
         FM_OUTPUT("is boundary known " << v.getId() << " " << p << endl);
      }
   }
   void setFar(const vertex &v)
   {
      status(&v) = FMStatus::F;
      scal Tv = L;
      setLs(v, Tv);
      tuple_data_bnd_t *pdbound = data_bnd.getData(v);
      if (pdbound) std::get<0>(*pdbound) = false;
   }
   bool isBoundaryLocalyUpdated(const vertex &v)
   {
      const tuple_data_bnd_t *pdbound = data_bnd.getData(v);
      if (pdbound)
      {
         if (std::get<0>(*pdbound))
            return true;
         else
            return false;
      }
      else
         return false;
   }
   bool isBoundary(const vertex &v)
   {
      const tuple_data_bnd_t *pdbound = data_bnd.getData(v);
      if (pdbound)
         return true;
      else
         return false;
   }

   template <typename ITER>
   void initBoundary(ITER b, ITER e)
   {
      tuple_data_bnd_t dbound;
      std::get<0>(dbound) = false;
      std::get<1>(dbound) = L;
      for (; b != e; ++b)
      {
         data_bnd.setData(**b) = dbound;
         FM_OUTPUT("is boundary " << (*b)->getId() << endl);
      }
   }

  private:
   const MESHINTERFACE &mi;
   DATAMANAGER<SCAL> &ls;
   std::function<SCAL(const vertex &)> Ffunc;
   const SCAL epsilon_ratio;
   data_bnd_t data_bnd;

   gls_t *pgls;
   vn_t *pvn;

   bool updatefp(const vertex &vx, const edge &e, const SCAL &epsilon)
   {
      SCAL valx;
      getGlobalLs(vx, valx);
      VECT px = getcoord(vx);
      SCAL Fx = Ffunc(vx);
      SCAL Trial;
      VECT GradTrial;
      VECT Grad;
      TRANSPORTED Trans;
      bool updated = false;
      const std::function<bool(const vertex *v)> is_known = [this](const vertex *v) {
         // return status(v);
         return (status(v) && status(v) != FMStatus::T);
      };
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
         tuple_data_bnd_t *pdbound = data_bnd.getData(vx);
         if (pdbound)
         {
            std::get<0>(*pdbound) = true;
         }
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

}  // namespace internal
}  // namespace xfastmarching
#endif
