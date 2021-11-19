/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <cassert>
#include <cmath>
#include <iostream>
// AOMD trellis
#include "AOMD_Internals.h"
#include "GEntity.h"
#include "mAOMD.h"
#include "mEntity.h"
#include "mVertex.h"
// xtensor
#include "xPoint.h"
#include "xTensor2.h"
#include "xVector.h"
// xmapping
#include "xMapping.h"
// xquadrature
#include "xIntegrator.h"
// xinterface
#include "xAOMDEntityUtil.h"
// xfem
#include "xDebug.h"
#include "xDomain.h"
#include "xGeomElem.h"
#include "xMappingBuilderHolder.h"

using namespace std;
using AOMD::mEntity;

namespace xfem
{
xGeomElem::xGeomElem(mEntity* e)
    : pent{e},
      ownedmapping{xMappingBuilderHolderSingleton::instance().buildMapping(*e)},
      mapping{ownedmapping.get()},
      ownedintegrator{new xquadrature::xGaussIntegrator(mapping->getType())},
      integrator{ownedintegrator.get()}
{
}

xGeomElem ::xGeomElem(mEntity* e, xmapping::xMapping* m, xquadrature::xIntegrator* i) : pent(e), mapping(m), integrator(i) {}

xGeomElem::xGeomElem(const xGeomElem& other)
    : pent{other.pent},
      Ori{other.Ori},
      order{other.order},
      uvw{other.uvw},
      CurrentIntegrationPoint{other.CurrentIntegrationPoint},
      NbIntegrationPoints{other.NbIntegrationPoints},
      Weight{other.Weight}
{
   if (other.ownedmapping)
   {
      ownedmapping.reset(xMappingBuilderHolderSingleton::instance().buildMapping(*pent));
      mapping = ownedmapping.get();
   }
   else
   {
      ownedmapping.reset();
      mapping = other.mapping;
   }
   if (other.ownedintegrator)
   {
      ownedintegrator.reset(new xquadrature::xGaussIntegrator(mapping->getType()));
      integrator = ownedintegrator.get();
   }
   else
   {
      ownedintegrator.reset();
      integrator = other.integrator;
   }
}

void swap(xGeomElem& l, xGeomElem& r) noexcept
{
   using std::swap;
   swap(l.pent, r.pent);
   swap(l.Ori, r.Ori);
   swap(l.ownedmapping, r.ownedmapping);
   swap(l.mapping, r.mapping);
   swap(l.ownedintegrator, r.ownedintegrator);
   swap(l.integrator, r.integrator);
   swap(l.order, r.order);
   swap(l.uvw, r.uvw);
   swap(l.CurrentIntegrationPoint, r.CurrentIntegrationPoint);
   swap(l.NbIntegrationPoints, r.NbIntegrationPoints);
   swap(l.Weight, r.Weight);
}

xGeomElem::xGeomElem(xGeomElem&& other) noexcept : xGeomElem() { swap(*this, other); }

xGeomElem& xGeomElem::operator=(xGeomElem other)
{
   swap(*this, other);
   return *this;
}

xtensor::xPoint xGeomElem ::getXYZ() const { return mapping->eval(uvw); }

xtensor::xPoint xGeomElem ::getCDGxyz() const { return mapping->eval(mapping->COG()); }

xtensor::xPoint xGeomElem ::getCDGuvw() const { return mapping->COG(); }

void xGeomElem ::setUVWForXYZ(const xtensor::xPoint& xyz)
{
   if (!mapping->invert(xyz, uvw))
   {
      std::cout << "Can't invert Mapping in " << __FILE__ << ", line : " << __LINE__ << std::endl;
      throw;
   }
}

void xGeomElem ::SetIntegrationPointNumberForDegree(const int& degtot)
{
   order = degtot;
   NbIntegrationPoints = integrator->nbIntegrationPoints(order);
}

int xGeomElem ::getZone() const
{
   int* pid = xDomain::get(*pent);
   return (pid) ? (*pid) : pent->getClassification()->tag();
}

void xGeomElem ::setUVW(const xtensor::xPoint& uvw_) { uvw = uvw_; }

void xGeomElem ::setUVW(int ipt)
{
   CurrentIntegrationPoint = ipt;
   integrator->iPoint(ipt, order, uvw(0), uvw(1), uvw(2), Weight);
}

void xGeomElem ::setUVWForVertex(int inod) { uvw = mapping->getMappingVertex(inod); }

xtensor::xVector<>& xGeomElem::PushBack(xtensor::xVector<>& inx) const
{
   mapping->pushBack(uvw, 1, &inx);
   return inx;
}

double xGeomElem::getMeasure()
{
   SetIntegrationPointNumberForDegree(1);
   double meas = 0.;
   const size_t nb = GetNbIntegrationPoints();
   for (size_t k = 0; k < nb; k++)
   {
      setUVW(int(k));
      meas += GetWeight() * GetDetJac();
   }
   return meas;
}

unsigned short iFace(const mEntity& e, const mEntity& border)
{
   auto l = e.getLevel();
   for (unsigned short ib = 0; ib < e.size(l - 1); ++ib)
      if (e.get(l - 1, ib) == &border) return ib;
   throw;
}

xtensor::xVector<> xGeomElem::normalVector(const mEntity& border) const
{
   return mapping->normalVector(iFace(*pent, border), uvw(0), uvw(1), uvw(2));
}

void xGeomElem::normalVector(mEntity* pborder, xtensor::xVector<>& n) const
{
   auto na = mapping->normalVector(iFace(*pent, *pborder), uvw(0), uvw(1), uvw(2));
   n = xtensor::xVector<>(na(0), na(1), na(2));
}

xgeom::xBoundingBox xGeomElem::getBoundingBox() const { return mapping->boundingBox(); }

void xGeomElem::GetBoundingBox(xtensor::xPoint& min, xtensor::xPoint& max) const
{
   xgeom::xBoundingBox bb = getBoundingBox();
   min = bb.min;
   max = bb.max;
}

bool xGeomElem::pointInElement(const xtensor::xPoint& xyz) const
{
   xtensor::xPoint uvwloc;
   return mapping->interiorCheck(xyz, uvwloc);
}

///[[deprecated]]
bool xGeomElem::PointInElement(const xtensor::xPoint& xyz) const { return pointInElement(xyz); }

double xGeomElem::getDetJac() const { return mapping->detJac(uvw) * Ori; }

double xGeomElem::GetDetJac() const { return getDetJac(); }

xtensor::xTensor2<>& xGeomElem::pushBackRight(xtensor::xTensor2<>& inx) const
{
   xtensor::xTensor2<> jInvx;
   mapping->jacInverse(uvw, jInvx);
   inx = inx * jInvx;
   return inx;
}

xtensor::xTensor2<>& xGeomElem::PushBackRight(xtensor::xTensor2<>& inx) const { return pushBackRight(inx); }

xtensor::xTensor2<>& xGeomElem::pushBackRightTranspose(xtensor::xTensor2<>& inx) const
{
   xtensor::xTensor2<> jInvx;
   mapping->jacInverse(uvw, jInvx);
   inx = inx * jInvx.transpose();
   return inx;
}

xtensor::xTensor2<>& xGeomElem::PushBackRightTranspose(xtensor::xTensor2<>& inx) const { return PushBackRightTranspose(inx); }

void xGeomElem::setIntegrator(xquadrature::xIntegrator* _integrator)
{
   ownedintegrator.reset();
   integrator = _integrator;
}

void xGeomElem::setMapping(xmapping::xMapping* m)
{
   ownedmapping.reset();
   mapping = m;
}

}  // namespace xfem
