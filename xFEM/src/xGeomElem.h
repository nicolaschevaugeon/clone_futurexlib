/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xGEOMETRICAL_ELEMENT_H
#define _xGEOMETRICAL_ELEMENT_H
#include <memory>
// aomd trellis
#include "mEntity.h"
// xquadrature
#include "xIntegrator.h"
// xtensor
#include "xPoint.h"
#include "xTensor2.h"
#include "xVector.h"

namespace xmapping
{
class xMapping;
}

namespace xfem
{
///  xGeomElem combines an entity, a xMapping, an xIntegrator and a current point it the reference configuration.
///! It is heavily used to evaluate functions inside elements either in reference space or in geometrical space at current point.
///! or push gradient of function with regard to localcoordinate to global coordinate via the mapping.
///! UVW refers to the current point in the reference configuration. By default uvw is set to {0.,0.,0.}
///! UVW can be set manually with setUVW(const xPoint &) or setUVWforXYZ(const xPoint &)
///! or it can be set using the xIntegrator by setting the ith point in the current integration rule setUVW(int)
class xGeomElem
{
  public:
   /// Construct a GeomElem from a mEntity.
   //! The mapping is constructed using   xMappingBuilderHolderSingleton::instance().buildMapping(*e)
   //! And the integrator is constructed by calling xquadrature::xGaussIntegrator.
   //! The mapping and the integrator are detroyed with the GeomElem or when calling setMapping
   //! and setIntegrator respectively
   xGeomElem(AOMD::mEntity* e);
   /// Construct a GeomElem, from a mapping, and an integrator
   //! xGeomElem don't take ownership of the mapping nor the integrator in this case.
   //! They must not be destroyed before they are reset or the xGeomElem is destroy.
   xGeomElem(AOMD::mEntity* e, xmapping::xMapping* m, xquadrature::xIntegrator* i);
   /// Copy constructor. Note that the mapping and the integrator are fully copied if they where owned by other.
   xGeomElem(const xGeomElem& other);
   xGeomElem(xGeomElem&& other) noexcept;
   xGeomElem& operator=(xGeomElem other);
   friend void swap(xGeomElem& l, xGeomElem& r) noexcept;
   /// return the pointer to the mesh entity.
   AOMD::mEntity* getEntity() const { return const_cast<AOMD::mEntity*>(pent); }
   /// return the local coordinates of the current point.
   const xtensor::xPoint& getUVW() const { return uvw; }
   /// return the global coordinates of the curreent point (using the mapping)
   xtensor::xPoint getXYZ() const;
   /// return the center of gravity of the ref element in local coordinate.
   xtensor::xPoint getCDGuvw() const;
   /// return the center of gravity of the ref element pushed to global coordinate.
   xtensor::xPoint getCDGxyz() const;
   /// computes volume (3D), area (2D) or length (1D), using the current integrator, with rule of order 0
   double getMeasure();
   /// return the exterior normal to the GeomElem, on the given boder.
   //! The border must one of the dim-1 element that bound the entitie's geom Elem or throw.
   xtensor::xVector<> normalVector(const AOMD::mEntity& pborder) const;
   /// [deprecated("Use xtensor::xVector<> normalVector(const mEntity& pborder) const"]]
   void normalVector(AOMD::mEntity* border, xtensor::xVector<>& n) const;
   /// return the bounding box of the mapping.
   xgeom::xBoundingBox getBoundingBox() const;
   /// [[deprecated("Use getBoundingBox() instead")]]
   void GetBoundingBox(xtensor::xPoint& min, xtensor::xPoint& max) const;
   /// return true if point xyz is inside the element (according to the mapping)
   bool pointInElement(const xtensor::xPoint& xyz) const;
   // [[deprecated("Use bool PointInElement(const xtensor::xPoint& xyz) const; instead")]]
   bool PointInElement(const xtensor::xPoint& xyz) const;
   double getDetJac() const;
   /// [[deprecated("Use getDetJac() instead"]]
   double GetDetJac() const;
   /// return the Weight of the current point, as set by the integration rule.
   const double& getWeight() const { return Weight; }
   /// [[deprecated("Use getWeight() instead"]]
   const double& GetWeight() const { return getWeight(); }
   unsigned int getNbIntegrationPoints() const { return NbIntegrationPoints; }
   /// [[deprecated("Use getNbIntegrationPoints() t instead"]]
   unsigned int GetNbIntegrationPoints() const { return getNbIntegrationPoints(); }
   /// get/set orientation of the GeomElem : 1 :  same as the entity, -1 :  reversed. the jacobian is multiply by the orientation
   /// when computing detJac. \note theses function seems quite useless and impose some book keeping. We probably should remove
   /// them.
   int getOrientation() const { return Ori; }
   void setOrientation(int s) { Ori = (s >= 0) ? 1 : -1; }
   /// [[deprecated("Use getOrientation instead"]]
   int GetOrientation() const { return Ori; }
   /// [[deprecated("Use setOrientation instead"]]
   void SetOrientation(int s) { Ori = (s >= 0) ? 1 : -1; }
   /// set the degree of the integrator.
   void setIntegrationPointNumberForDegree(const int& degtot);
   ///  [[deprecated("Use setIntegrationPointNumberForDegree(const int& degtot) instead")]]
   void SetIntegrationPointNumberForDegree(const int& degtot);
   /// return the order of the integration rule (as set by ) setIntegrationPointNumberForDegree
   unsigned int getOrder() { return order; }
   /// set the current local point to the ipt point of the integration rule.
   //! ipt must be lower than the number of points in the current integration rule (getNbIntegrationPoints());
   void setUVW(int ipt);
   /// set the current point local coordinates
   void setUVW(const xtensor::xPoint& uvw_);
   /// set the current point by inverting the mapping at xyz. Can fail if xyz is not on the element.
   void setUVWForXYZ(const xtensor::xPoint& xyz);
   /// set UVW at the mapping vertex inod (see xMapping)
   void setUVWForVertex(int inod);
   /// return the current point number in the integration rule.
   int getCurrentIntegrationPointId() const { return CurrentIntegrationPoint; }
   xtensor::xVector<>& pushBack(xtensor::xVector<>& inx) const;
   xtensor::xTensor2<>& pushBackRight(xtensor::xTensor2<>& inx) const;
   xtensor::xTensor2<>& pushBackRightTranspose(xtensor::xTensor2<>& inx) const;
   /// [[deprecated("Use pushBack(xtensor::xVector<>& inx) instead")]]
   xtensor::xVector<>& PushBack(xtensor::xVector<>& inx) const;
   /// [[deprecated("Use pushBackRight(xtensor::xVector<>& inx) instead")]]
   xtensor::xTensor2<>& PushBackRight(xtensor::xTensor2<>& inx) const;
   /// [[deprecated("Use pushBackRightTranspose(xtensor::xVector<>& inx) instead")]]
   xtensor::xTensor2<>& PushBackRightTranspose(xtensor::xTensor2<>& inx) const;
   ///  return the element's zone tag.
   int getZone() const;
   /// return the mapping
   const xmapping::xMapping* getMapping() const { return mapping; }
   /// set the mapping. if mapping set this way, the user is responsible for the pointed to lifetime.
   void setMapping(xmapping::xMapping* m);
   /// return the quadrature.
   const xquadrature::xIntegrator* getIntegrator() const { return integrator; }
   /// set the integrator. if setted this way, the user is responsible for the life time.
   void setIntegrator(xquadrature::xIntegrator* integrator);

  private:
   // \note : need to be a pointer to cope with operator =.
   const AOMD::mEntity* pent = nullptr;
   int Ori = 1;
   // \note the mapping and the integrator can be owned or not by the xGeomElem.
   // so we have a unique pointer ownedmapping set to the mapping if xGeomElem owns it,
   // and mapping point to the ownedmapping target in this case.
   // if xGeomElem don't own it,  ownedmapping is set to nullptr and mapping point to what ever the user want,
   // but wat is pointed to  won't be destroyed at detruction of xGeomElem.
   // the same idiom is used for xintegrator.
   std::unique_ptr<xmapping::xMapping> ownedmapping = nullptr;
   xmapping::xMapping* mapping = nullptr;
   std::unique_ptr<xquadrature::xIntegrator> ownedintegrator = nullptr;
   xquadrature::xIntegrator* integrator = nullptr;
   int order = 0;
   xtensor::xPoint uvw;
   int CurrentIntegrationPoint = 0;
   unsigned int NbIntegrationPoints = 1;
   double Weight = 1.;
   xGeomElem() = default;
};

}  // namespace xfem

#endif
