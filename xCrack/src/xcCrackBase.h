/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef _XCCRACKBASE_
#define _XCCRACKBASE_
#include <functional>
#include <vector>

#include "xEntityFilter.h"
#include "xTensor2.h"
namespace AOMD
{
class mEntity;
}
namespace xfem
{
class xGeomElem;
class xLevelSet;
}  // namespace xfem

namespace xcrack
{
/// In order to update and later replace the old lCrack class, A Virtual base class xcCrackBase is build in the following.
/*! The xcCrackBase contain the virtual definitions of all the functions necessary to define enrichment function.!*/
class xcCrackBase
{
  public:
   virtual ~xcCrackBase() = default;
   /// A filter that return true if the entities in entry has its support completly cut by the crack.
   /*! Primarly used for Heaviside enrichment
   !*/
   virtual bool supportsCutthruByCrack(AOMD::mEntity* e) const = 0;
   virtual xfem::xEntityFilter supportsCutthruByCrackFilter() const
   {
      return std::bind1st(std::mem_fun(&xcrack::xcCrackBase::supportsCutthruByCrack), this);
   }
   /// This evaluator is meant to decide if an entity is above or below the crack surface.
   /*! Return 1 on top of the crack surface (lsn > 0), zero or -1 otherwise.
       It uses the usual xEval convention, geo_appro store the entity on which the approximation is defined, and geo_integ store
   the entity on which the information is evaluated. The result can be overriden to either -1 or 1  if member overrideside =-1 or
   1
   !*/
   virtual int sideOf(const xfem::xGeomElem* geo_appro, const xfem::xGeomElem* geo_integ) const = 0;
   /// Get distance from  the lip and angle;
   /* From  emesh entity e and coordinate in the reference frame on this element, return r = distance to the crack lip and
      theta direction angle, 0 behing the crack plan
   */
   virtual void getLocalCoords(AOMD::mEntity* e, xtensor::xPoint& uvw, double& r, double& theta) const = 0;

   virtual std::vector<double> getLstValues(AOMD::mEntity* e) const = 0;
   virtual std::vector<double> getLsnValues(AOMD::mEntity* e) const = 0;

   virtual void getCrackAxis(AOMD::mEntity* e, xtensor::xTensor2<>& base) const = 0;
   virtual void getCrackOrthoAxis(AOMD::mEntity* e, xtensor::xTensor2<>& base) const = 0;

   /// Transform the coordinate of a vector expressed in the crack coordinate system (t, n, t^n) to the global coordinate system
   /*!
     Entity e and uvw (coordinate of the point on which the vector is computed in the reference element of e)
     are used to compute the frame change.
     !*/
   virtual void localToGlobal(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& local) const = 0;
   /// From one element on which lsn and lst must be defined, return the crack tip bases and their gradient.
   virtual void getLocalCurv(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xVector<>& e1, xtensor::xVector<>& e2,
                             xtensor::xVector<>& e3, xtensor::xTensor2<>& curv1, xtensor::xTensor2<>& curv2,
                             xtensor::xTensor2<>& curv3) const = 0;
   /// Transform the coordinate of a tensor2 expressed in the crack coordinate system (t, n, t^n) to the global coordinate system
   /*!
     Entity e and uvw (coordinate of the point on which the Tensor  is computed in the reference element of e)
     are used to compute the frame change. The "smoothed gradient" is use here
     !*/
   virtual void localToGlobal(AOMD::mEntity* e, const xtensor::xPoint& uvw, xtensor::xTensor2<>& local) const = 0;

   virtual xfem::xMesh* getMeshCrackFront() = 0;
   virtual const xfem::xMesh& getMeshCrackFront() const = 0;
   virtual xfem::xLevelSet const* getFieldt() const = 0;
   virtual xfem::xLevelSet const* getFieldn() const = 0;
   virtual xfem::xMesh* getMesh() = 0;
   virtual const xfem::xMesh& getMesh() const = 0;
   virtual double getFrontDistance(AOMD::mVertex* v) const = 0;
   std::string label;
   int overrideside = 0;
   virtual void getAsymptoticFields(const xfem::xGeomElem* geo_appro, int mode, xtensor::xTensor2<>& stress_aux,
                                    xtensor::xTensor2<>& grad_disp_aux) const = 0;
};

}  // namespace xcrack
#endif
