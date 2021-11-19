/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _OPERATORS_H
#define _OPERATORS_H

#include <vector>

#include "xApproxFctPtr.h"
#include "xApproxFunction.h"
#include "xField.h"
#include "xGeomElem.h"

namespace xfem
{
template <class UnaryOperator>
class xValOperator
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   typedef typename UnaryOperator::argument_type argument_type;
   xValOperator() : funct() {}
   const UnaryOperator funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      vals.reserve(vals.size() + f->size());
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);  // FF are assumed double valued
   }
};

template <class UnaryOperator>
class xGradOperator
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   const UnaryOperator funct;
   xGradOperator() : funct() {}
   xGradOperator(UnaryOperator& op) : funct(op) {}
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      vals.reserve(vals.size() + f->size());
      xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);  // FF are assumed double valued
   }
};

template <class UnaryOperator>
class xGradOperatorAxisym
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   const UnaryOperator funct;
   xGradOperatorAxisym() : funct() {}
   xGradOperatorAxisym(UnaryOperator& op) : funct(op) {}
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      vals.reserve(vals.size() + f->size());
      xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

      std::vector<xtensor::xVector<>> vals_tmp;
      vals_tmp.reserve(vals.size() + f->size());
      xField<>::getFF(f->begin(), f->end(), vals_tmp, geo_appro, geo_integ, xtool::xIdentity<xtensor::xVector<>>());

      auto uvw = geo_integ->getXYZ();
      double r = uvw(0);
      for (int i = 0; i < vals.size() / 2; ++i)
      {
         vals[i](2, 2) = vals_tmp[i](0) / r;
      }
   }
};

template <class UnaryOperator>
class xGradLocalOperator
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   const UnaryOperator funct;
   xGradLocalOperator() : funct() {}
   xGradLocalOperator(UnaryOperator& op) : funct(op) {}
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      xField<>::getGradLocalFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);
   }
};

template <class UnaryOperator>
class xGradLocalOperatorAxisym
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   const UnaryOperator funct;
   xGradLocalOperatorAxisym() : funct() {}
   xGradLocalOperatorAxisym(UnaryOperator& op) : funct(op) {}
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      xField<>::getGradLocalFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

      std::vector<xtensor::xVector<>> vals_tmp;
      vals_tmp.reserve(vals.size() + f->size());
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, xtool::xIdentity<xtensor::xVector<>>());

      auto uvw = geo_integ->getXYZ();
      double r = uvw(0);
      for (int i = 0; i < vals.size() / 2; ++i)
      {
         vals[i](2, 2) = vals_tmp[i](0) / r;
      }
   }
};

class xProjectOnVector : public std::unary_function<xtensor::xVector<>, double>
{
  public:
   xProjectOnVector(const xEval<xtensor::xVector<>>& _eval_vec) : eval_vec(_eval_vec){};
   void setGeom(xGeomElem* geo_appro, xGeomElem* geo_integ) { eval_vec(geo_appro, geo_integ, vec); }

   double operator()(xtensor::xVector<>& in) const { return in(0) * vec(0) + in(1) * vec(1) + in(2) * vec(2); }

  private:
   const xEval<xtensor::xVector<>>& eval_vec;
   xtensor::xVector<> vec;
};

class xProjectOnPlane : public std::unary_function<xtensor::xVector<>, xtensor::xVector<>>
{
  public:
   xProjectOnPlane(xEval<xtensor::xVector<>>& _eval_vec) : eval_vec(_eval_vec){};
   void setGeom(xGeomElem* geo_appro, xGeomElem* geo_integ) { eval_vec(geo_appro, geo_integ, vec); }

   xtensor::xVector<> operator()(xtensor::xVector<>& in) const
   {
      double a = in * vec;
      return in - vec * a;
   }

  private:
   const xEval<xtensor::xVector<>>& eval_vec;
   xtensor::xVector<> vec;
};

template <class UnaryOperator>
class xValOperatorWithParam
{
  public:
   typedef typename UnaryOperator::result_type result_type;

  public:
   xValOperatorWithParam(UnaryOperator& f) : funct(f) {}
   UnaryOperator& funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      funct.setGeom(geo_appro, geo_integ);
      vals.reserve(vals.size() + f->size());
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);
   }
};

template <class UnaryOperator>
class xGradOperatorWithParam
{
  public:
   typedef typename UnaryOperator::result_type result_type;

  public:
   xGradOperatorWithParam(UnaryOperator& f) : funct(f) {}
   UnaryOperator& funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      funct.setGeom(geo_appro, geo_integ);
      vals.reserve(vals.size() + f->size());
      xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);
   }
};

/// \brief To compute the jump of the approximation across geo_integ
template <class UnaryOperator>
class xValJumpOperator
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   xValJumpOperator(const char* side_tag_name = "side_tag")
       : funct(), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name))
   {
   }
   const UnaryOperator funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      const int rs = vals.size() + f->size();
      std::vector<result_type> valsNeg = vals;
      vals.reserve(rs);
      valsNeg.reserve(rs);

      // Positive side
      geo_appro->getEntity()->attachInt(side_tag, 1);
      geo_integ->getEntity()->attachInt(side_tag, 1);
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

      // Negative side
      geo_appro->getEntity()->attachInt(side_tag, -1);
      geo_integ->getEntity()->attachInt(side_tag, -1);
      xField<>::getFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct, true);

      // Untag side:
      geo_appro->getEntity()->deleteData(side_tag);
      geo_integ->getEntity()->deleteData(side_tag);

      // Reset Storage before leaving (test)
      xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

      // Compute jump :
      // std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::minus<result_type>()); //Bug ?
      for (unsigned int i = 0; i < vals.size(); ++i)
      {
         vals[i] -= valsNeg[i];
      }
   }

  private:
   unsigned int side_tag;
};

template <class UnaryOperator>
class xValJumpOperatorWithParam
{
  public:
   typedef typename UnaryOperator::result_type result_type;
   xValJumpOperatorWithParam(UnaryOperator& f, const char* side_tag_name = "side_tag")
       : funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name))
   {
   }
   UnaryOperator& funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      const int rs = vals.size() + f->size();
      std::vector<result_type> valsNeg = vals;
      vals.reserve(rs);
      valsNeg.reserve(rs);

      // Positive side
      geo_appro->getEntity()->attachInt(side_tag, 1);
      geo_integ->getEntity()->attachInt(side_tag, 1);
      funct.setGeom(geo_appro, geo_integ);
      xField<>::getFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

      // Negative side
      geo_appro->getEntity()->attachInt(side_tag, -1);
      geo_integ->getEntity()->attachInt(side_tag, -1);
      funct.setGeom(geo_appro, geo_integ);
      xField<>::getFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct, true);

      // Untag side:
      geo_appro->getEntity()->deleteData(side_tag);
      geo_integ->getEntity()->deleteData(side_tag);

      // Reset Storage before leaving (test)
      xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

      // Compute jump :
      // std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::minus<result_type>()); //Bug ?
      for (unsigned int i = 0; i < vals.size(); ++i)
      {
         vals[i] -= valsNeg[i];
      }
   }

  private:
   unsigned int side_tag;
};

/// \brief To compute the mean of the gradient of the approximation across geo_integ
/// \todo Faire la fin de la moyenne avec un transform... je ne comprends pas pq ca ne fct pas !
template <class UnaryOperator>
class xGradMeanOperatorWithParam
{
  public:
   typedef typename UnaryOperator::result_type result_type;

  public:
   xGradMeanOperatorWithParam(UnaryOperator& f, const char* side_tag_name = "side_tag")
       : funct(f), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name))
   {
   }
   UnaryOperator& funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      const int rs = vals.size() + f->size();
      std::vector<result_type> valsNeg = vals;
      vals.reserve(rs);
      valsNeg.reserve(rs);

      // Positive side
      geo_appro->getEntity()->attachInt(side_tag, 1);
      geo_integ->getEntity()->attachInt(side_tag, 1);
      funct.setGeom(geo_appro, geo_integ);
      xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

      // Negative side
      geo_appro->getEntity()->attachInt(side_tag, -1);
      geo_integ->getEntity()->attachInt(side_tag, -1);
      funct.setGeom(geo_appro, geo_integ);
      xField<>::getGradFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct, true);

      // Untag side:
      geo_appro->getEntity()->deleteData(side_tag);
      geo_integ->getEntity()->deleteData(side_tag);

      // Reset Storage before leaving (test)
      xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

      // Compute mean :
      // std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::plus<result_type>());

      // A faire fonctionner !!!! (et virer la boucle for) : JE NE COMPRENDS PAS CE QUI NE VA PAS DANS LE TRANSFORM....
      //     std::transform(vals.begin(),vals.end(),vals.begin(),bind2nd(std::multiplies<result_type>(),0.5));
      // for(int i =0; i<vals.size();++i) vals[i]*=0.5;

      for (unsigned int i = 0; i < vals.size(); ++i) vals[i] = (vals[i] + valsNeg[i]) * 0.5;
   }

  private:
   unsigned int side_tag;
};

template <class UnaryOperator>
class xGradMeanOperator
{
  public:
   typedef typename UnaryOperator::result_type result_type;

  public:
   xGradMeanOperator(const char* side_tag_name = "side_tag")
       : funct(), side_tag(AOMD::AOMD_Util::Instance()->lookupMeshDataId(side_tag_name))
   {
   }
   const UnaryOperator funct;
   void eval(femFcts_t* f, xGeomElem* geo_appro, xGeomElem* geo_integ, std::vector<result_type>& vals) const
   {
      const int rs = vals.size() + f->size();
      std::vector<result_type> valsNeg = vals;
      vals.reserve(rs);
      valsNeg.reserve(rs);

      // Positive side
      geo_appro->getEntity()->attachInt(side_tag, 1);
      geo_integ->getEntity()->attachInt(side_tag, 1);
      xField<>::getGradFF(f->begin(), f->end(), vals, geo_appro, geo_integ, funct);

      // Negative side
      geo_appro->getEntity()->attachInt(side_tag, -1);
      geo_integ->getEntity()->attachInt(side_tag, -1);
      xField<>::getGradFF(f->begin(), f->end(), valsNeg, geo_appro, geo_integ, funct, true);

      // Untag side:
      geo_appro->getEntity()->deleteData(side_tag);
      geo_integ->getEntity()->deleteData(side_tag);

      // Reset Storage before leaving (test)
      xField<>::resetStorage(f->begin(), f->end(), geo_appro, geo_integ);

      // Compute mean :
      // std::transform(vals.begin(),vals.end(),valsNeg.begin(),vals.begin(),std::plus<result_type>());

      // A faire fonctionner !!!! (et virer la boucle for) : JE NE COMPRENDS PAS CE QUI NE  VA PAS DANS LE TRANSFORM....
      // for(int i =0; i<vals.size();++i) vals[i]*=0.5;

      for (int i = 0; i < vals.size(); ++i) vals[i] = (vals[i] + valsNeg[i]) * 0.5;
   }

  private:
   unsigned int side_tag;
};

// If needed We still need to code the HessianOperator
}  // namespace xfem

#endif
