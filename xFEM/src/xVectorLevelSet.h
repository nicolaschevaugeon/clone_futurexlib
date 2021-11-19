/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#ifndef XVECTORLEVELSET
#define XVECTORLEVELSET

#include <set>
#include <string>

// xinterface aomd
#include "xAttachedDataManagerAOMD.h"
// xfem
#include "xElement.h"
#include "xEntityFilter.h"
#include "xGetSupport.h"
#include "xLevelSet.h"
#include "xPartition.h"
#include "xVector.h"

namespace xfem
{
class xVectorLevelSetData
{
  public:
   xVectorLevelSetData() = delete;
   xVectorLevelSetData(const int& _inout, const xtensor::xVector<>& _directiontoiso, const bool exist_ = true)
       : inout(_inout), directiontoiso(_directiontoiso), exist(exist_)
   {
   }
   int getInOut() const { return inout; }
   bool doExist() const { return exist; }
   xtensor::xVector<> getDirectionToIso() const { return directiontoiso; }

  private:
   // \note I unconstify to allow for default constructor, needed by DataManager. Could put back const and remove default
   // constructor
   //  if DataManager have the possibility to inplace construct....
   const int inout;
   const xtensor::xVector<> directiontoiso;
   const bool exist;
};

inline bool compLessThanVLSDataInOut(const xVectorLevelSetData& left, const xVectorLevelSetData& right)
{
   return (left.getInOut() < right.getInOut());
}

inline bool compLessThanVLSDataInOut(const xVectorLevelSetData* left, const xVectorLevelSetData* right)
{
   if (left && right)
   {
      return compLessThanVLSDataInOut(*left, *right);
   }
   else
      throw;
}

typedef std::function<xVectorLevelSetData(const xtensor::xPoint&)> xPointToVectorLevelSetData;

class xAnalyticalVectorLevelSetCylinder
{
  public:
   xAnalyticalVectorLevelSetCylinder(const xtensor::xPoint& p, xtensor::xVector<> d, const double& r);
   xVectorLevelSetData operator()(const xtensor::xPoint&) const;

  private:
   const xtensor::xPoint point_on_axis;
   const xtensor::xVector<> axis;
   const double radius;
   xtensor::xVector<> ortho;
};

class xAnalyticalVectorLevelSetSphere
{
  public:
   xAnalyticalVectorLevelSetSphere(const xtensor::xPoint& center, const double& radius);
   xVectorLevelSetData operator()(const xtensor::xPoint&) const;

  private:
   const xtensor::xPoint center;
   const double radius;
};

class xAnalyticalVectorLevelSetGrownSegment
{
  public:
   xAnalyticalVectorLevelSetGrownSegment(const xtensor::xPoint& p1_, const xtensor::xPoint& p2_, const double& growth_);
   xVectorLevelSetData operator()(const xtensor::xPoint& p) const;

  private:
   const xtensor::xPoint p1, p2;
   const xAnalyticalVectorLevelSetSphere sphere1;
   const xAnalyticalVectorLevelSetSphere sphere2;
   const xAnalyticalVectorLevelSetCylinder cyl;
   xtensor::xVector<> dir;
   double length;
};

template <class AVLS0, class AVLS1>
class xAnalyticalVectorLevelSetUnion
{
  public:
   xAnalyticalVectorLevelSetUnion(const AVLS0& vls0, const AVLS1& vls1);
   xVectorLevelSetData operator()(const xtensor::xPoint& p) const;

  private:
   AVLS0 vls0;
   AVLS1 vls1;
};

template <class AVLS0, class AVLS1>
xAnalyticalVectorLevelSetUnion<AVLS0, AVLS1>::xAnalyticalVectorLevelSetUnion(const AVLS0& vls0_, const AVLS1& vls1_)
    : vls0{vls0_}, vls1{vls1_}
{
}

template <class AVLS0, class AVLS1>
xVectorLevelSetData xAnalyticalVectorLevelSetUnion<AVLS0, AVLS1>::operator()(const xtensor::xPoint& p) const
{
   const xVectorLevelSetData vls0p = vls0(p);
   const xVectorLevelSetData vls1p = vls1(p);
   const int inout0 = vls0p.getInOut();
   const int inout1 = vls1p.getInOut();
   xtensor::xVector<double> d0(vls0p.getDirectionToIso());
   xtensor::xVector<double> d1(vls1p.getDirectionToIso());
   if (inout0 >= 0 && inout1 >= 0)
   {
      return (d0.mag() < d1.mag()) ? vls0p : vls1p;
   }
   if (inout0 < 0 && inout1 >= 0) return vls0p;
   if (inout1 < 0 && inout0 >= 0) return vls1p;
   const xtensor::xPoint p0 = p + d0;
   const xVectorLevelSetData vls1p0 = vls1(p0);
   if (vls1p0.getInOut() >= 0) return xVectorLevelSetData(-1, d0);
   return xVectorLevelSetData(-1, d0 + vls1p0.getDirectionToIso());
}

class xVectorLevelSet
{
  public:
   template <typename VERTEXITERATOR>
   xVectorLevelSet(const VERTEXITERATOR& _vbegin, const VERTEXITERATOR& _vend, xPointToVectorLevelSetData vlseval);
   template <typename VERTEXITERATOR, typename GETCLOSEST>
   xVectorLevelSet(const VERTEXITERATOR& _vbegin, const VERTEXITERATOR& _vend, const xLevelSet& ls, GETCLOSEST& getclosestpoint,
                   const double& shift);
   ~xVectorLevelSet();
   const xVectorLevelSetData* operator()(const AOMD::mVertex& v) const;
   xtensor::xVector<> getValVector(AOMD::mEntity* e, const xtensor::xPoint& uvw) const;
   double getValInOut(AOMD::mEntity* e, const xtensor::xPoint& uvw) const;
   // void exportGmsh(const std::string &filename) const;
  private:
   xinterface::aomd::xAttachedDataManagerAOMD<xVectorLevelSetData> data;
};

template <class UnaryOperator>
class xEvalVectorLevelSetVector : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalVectorLevelSetVector(const xVectorLevelSet& vls_) : vls(vls_) {}
   xEvalVectorLevelSetVector(const xVectorLevelSet& vls_, const UnaryOperator& _funct) : vls(vls_), funct(_funct) {}
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v;
      v = vls.getValVector(appro->getEntity(), appro->getUVW());
      result = funct(v);
   }

  private:
   const xVectorLevelSet& vls;
   UnaryOperator funct;
};
template <class UnaryOperator>
class xEvalVectorLevelSetInOut : public xEval<typename UnaryOperator::result_type>
{
  public:
   typedef typename xEval<typename UnaryOperator::result_type>::result_type result_type;

   xEvalVectorLevelSetInOut(const xVectorLevelSet& vls_) : vls(vls_) {}
   xEvalVectorLevelSetInOut(const xVectorLevelSet& vls_, const UnaryOperator& _funct) : vls(vls_), funct(_funct) {}
   void operator()(const xGeomElem* appro, const xGeomElem* integ, result_type& result) const override
   {
      typename UnaryOperator::argument_type v;
      v = vls.getValInOut(appro->getEntity(), appro->getUVW());
      result = funct(v);
   }

  private:
   const xVectorLevelSet& vls;
   UnaryOperator funct;
};

}  // namespace xfem

#include "xVectorLevelSet_imp.h"

#endif
