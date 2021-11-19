/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <numeric>

// aomd
#include "mEntity.h"
#include "mFace.h"
#include "mHex.h"
#include "mVertex.h"
// xfem
#include "xApproxFunction.h"
#include "xDebug.h"
#include "xElement.h"
// xmapping
#include "xLagrangeMapping.h"
#include "xReferenceElement.h"

namespace xfem
{
using AOMD::mEntity;
using AOMD::mVertex;
using std::cerr;

xElement::xElement(const mEntity* e) : dim(e->getLevel()), type(e->getType())
{
   assert(e->getType() == mEntity::VERTEX || e->getType() == mEntity::EDGE || e->getType() == mEntity::TET ||
          e->getType() == mEntity::TRI || e->getType() == mEntity::QUAD || e->getType() == mEntity::HEX);

   for (int j = 0; j < e->size(0); j++)
   {
      mVertex* v = static_cast<mVertex*>(e->get(0, j));
      p.push_back(v->point());
   }
   switch (dim)
   {
      case 3:
      {
         // for the jacobian matrix, the first index is uvw
         // the second index                         is xyz
         jac(0, 0) = p[1](0) - p[0](0);
         jac(0, 1) = p[1](1) - p[0](1);
         jac(0, 2) = p[1](2) - p[0](2);
         jac(1, 0) = p[2](0) - p[0](0);
         jac(1, 1) = p[2](1) - p[0](1);
         jac(1, 2) = p[2](2) - p[0](2);
         jac(2, 0) = p[3](0) - p[0](0);
         jac(2, 1) = p[3](1) - p[0](1);
         jac(2, 2) = p[3](2) - p[0](2);
         detjac = jac(0, 0) * (jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1)) +
                  jac(1, 0) * (jac(2, 1) * jac(0, 2) - jac(0, 1) * jac(2, 2)) +
                  jac(2, 0) * (jac(0, 1) * jac(1, 2) - jac(1, 1) * jac(0, 2));
      }
      break;
      case 2:
      {
         if (type == mEntity::TRI)
         {
            jac(0, 0) = p[1](0) - p[0](0);
            jac(0, 1) = p[1](1) - p[0](1);
            jac(0, 2) = p[1](2) - p[0](2);
            jac(1, 0) = p[2](0) - p[0](0);
            jac(1, 1) = p[2](1) - p[0](1);
            jac(1, 2) = p[2](2) - p[0](2);
            double d3 = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
            double d2 = jac(0, 2) * jac(1, 0) - jac(0, 0) * jac(1, 2);
            double d1 = jac(0, 1) * jac(1, 2) - jac(0, 2) * jac(1, 1);
            detjac = std::sqrt(d1 * d1 + d2 * d2 + d3 * d3);
            jac(2, 0) = d1 / detjac;
            jac(2, 1) = d2 / detjac;
            jac(2, 2) = d3 / detjac;
         }
         else
         { /*QUAD*/
            // Do nothing : we do not want to implement things such as propagation.
         }
      }
      break;
      case 1:
      {
         jac(0, 0) = 0.5 * (p[1](0) - p[0](0));
         jac(0, 1) = 0.5 * (p[1](1) - p[0](1));
         jac(0, 2) = 0.5 * (p[1](2) - p[0](2));
         detjac = sqrt(jac(0, 0) * jac(0, 0) + jac(0, 1) * jac(0, 1) + jac(0, 2) * jac(0, 2));
         int vrai = 0;
         for (int i = 0; i < 3; i++)
         {
            if (jac(0, i) == 0)
            {  // if this component of the normal vector is zero
               vrai = 1;
               for (int j = 0; j < 3; j++)
               {
                  if (j == i)
                     jac(1, j) = 1;  // then the component of the second normal must be one,
                  else
                     jac(1, j) = 0;  // and the other components of the second normal are zero.
               }
               continue;
            }
         }
         if (!vrai)
         {  // looking for the second normal vector in the plane z = 0
            double temp = sqrt(jac(0, 0) * jac(0, 0) + jac(0, 1) * jac(0, 1));
            jac(1, 0) = -jac(0, 1) / temp;
            jac(1, 1) = jac(0, 0) / temp;
            jac(1, 2) = 0;
         }
         // The third normal vector
         jac(2, 0) = (jac(0, 1) * jac(1, 2) - jac(0, 2) * jac(1, 1)) / detjac;
         jac(2, 1) = (jac(0, 2) * jac(1, 0) - jac(0, 0) * jac(1, 2)) / detjac;
         jac(2, 2) = (jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0)) / detjac;
      }
      break;
      case 0:
      {
         jac(0, 0) = jac(1, 1) = jac(2, 2) = 1.;
         jac(0, 1) = jac(0, 1) = jac(0, 2) = jac(0, 2) = jac(1, 2) = jac(2, 1) = 0.;
         detjac = 1.;
      }
      break;
      default:
         cerr << "valid only for EDGE, TRI and TET\n";
         assert(0);
         break;
   }

   invjac(0, 0) = (jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1)) / detjac;
   invjac(1, 0) = (jac(1, 2) * jac(2, 0) - jac(1, 0) * jac(2, 2)) / detjac;
   invjac(2, 0) = (jac(1, 0) * jac(2, 1) - jac(2, 0) * jac(1, 1)) / detjac;

   invjac(0, 1) = (jac(2, 1) * jac(0, 2) - jac(0, 1) * jac(2, 2)) / detjac;
   invjac(1, 1) = (jac(0, 0) * jac(2, 2) - jac(2, 0) * jac(0, 2)) / detjac;
   invjac(2, 1) = (jac(2, 0) * jac(0, 1) - jac(0, 0) * jac(2, 1)) / detjac;

   invjac(0, 2) = (jac(0, 1) * jac(1, 2) - jac(1, 1) * jac(0, 2)) / detjac;
   invjac(1, 2) = (jac(1, 0) * jac(0, 2) - jac(0, 0) * jac(1, 2)) / detjac;
   invjac(2, 2) = (jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1)) / detjac;

   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         tinvjac(i, j) = invjac(j, i);
      }
   }
}

void xElement::getFFTri(const double& u, const double& v, std::vector<double>& ff)
{
   ff.push_back(1. - u - v);
   ff.push_back(u);
   ff.push_back(v);
}
void xElement::getFFQuadSimplex(const double& u, const double& v, std::vector<double>& ff)
{
   ff.push_back(SimplexApproxFunctionQuadrilateral(1, u, v));
   ff.push_back(SimplexApproxFunctionQuadrilateral(2, u, v));
   ff.push_back(SimplexApproxFunctionQuadrilateral(3, u, v));
   ff.push_back(SimplexApproxFunctionQuadrilateral(4, u, v));
}
void xElement::getFFHexSimplex(const double& u, const double& v, const double& w, std::vector<double>& ff)
{
   ff.push_back(SimplexApproxFunctionHexahedron(1, u, v, w));
   ff.push_back(SimplexApproxFunctionHexahedron(2, u, v, w));
   ff.push_back(SimplexApproxFunctionHexahedron(3, u, v, w));
   ff.push_back(SimplexApproxFunctionHexahedron(4, u, v, w));
   ff.push_back(SimplexApproxFunctionHexahedron(5, u, v, w));
   ff.push_back(SimplexApproxFunctionHexahedron(6, u, v, w));
   ff.push_back(SimplexApproxFunctionHexahedron(7, u, v, w));
   ff.push_back(SimplexApproxFunctionHexahedron(8, u, v, w));
}
void xElement::getFFEdge(const double& u, std::vector<double>& ff)
{
   ff.push_back(0.5 * (1. - u));
   ff.push_back(0.5 * (1. + u));
}
void xElement::getFFTet(const double& u, const double& v, const double& w, std::vector<double>& ff)
{
   ff.push_back(1. - u - v - w);
   ff.push_back(u);
   ff.push_back(v);
   ff.push_back(w);
}

void xElement::getFF(std::vector<double>& ff)
{
   //  ... cette fonction definit les fonctions de formes a partir des coordonnees locales
   ff.reserve(dim + 1);
   if (dim == 3)
   {
      if (type == mEntity::TET)
         getFFTet(uvw(0), uvw(1), uvw(2), ff);
      else
      {
         ff.reserve(8);
         getFFHexSimplex(uvw(0), uvw(1), uvw(2), ff);
      }
   }
   else if (dim == 2)
   {
      if (type == mEntity::TRI)
         getFFTri(uvw(0), uvw(1), ff);
      else
      {
         ff.reserve(4);
         getFFQuadSimplex(uvw(0), uvw(1), ff);
      }
   }
   else
      getFFEdge(uvw(0), ff);
}

// ICI, on suppose qu on connait les valeurs des fonctions de forme au point considere
double xElement::getInterpoSca(const std::vector<double>& scas)
{
   //  ... calcul d un champ scalaire interpole avec les fonctions de formes
   std::vector<double> ff;
   getFF(ff);
   return std::inner_product(ff.begin(), ff.end(), scas.begin(), 0.0);
}

// ICI, on suppose qu on connait les valeurs des fonctions de forme au point considere
xtensor::xVector<> xElement::getInterpoVec(const std::vector<xtensor::xVector<>>& vecs)
{
   //  ... calcul d un champ vectoriel interpole avec les fonctions des formes
   std::vector<double> ff;
   getFF(ff);
   xtensor::xVector<> out(std::inner_product(vecs.begin(), vecs.end(), ff.begin(), xtensor::xVector<>(0., 0., 0)));
   return out;
}

void xElement::getGradFF(std::vector<xtensor::xVector<>>& gradff)
{
   const bool debug = xdebug_flag;
   //  ... cette fonction definit les gradients des fonctions de formes dans la base globale
   gradff.clear();
   gradff.reserve(dim + 1);
   if (dim == 3)
   {
      if (type == mEntity::TET)
      {
         gradff.push_back(invjac * xtensor::xVector<>(-1., -1., -1.));
         gradff.push_back(invjac * xtensor::xVector<>(1., 0., 0.));
         gradff.push_back(invjac * xtensor::xVector<>(0., 1., 0.));
         gradff.push_back(invjac * xtensor::xVector<>(0., 0., 1.));
      }
      else
         assert(0); /*Not implemented... kicked out!*/
   }
   else if (dim == 2)
   {
      if (type == mEntity::TRI)
      {
         gradff.push_back(invjac * xtensor::xVector<>(-1., -1., 0.));
         gradff.push_back(invjac * xtensor::xVector<>(1., 0., 0.));
         gradff.push_back(invjac * xtensor::xVector<>(0., 1., 0.));
      }
      else
         assert(0); /*Not implemented... kicked out!*/
   }
   else
   {
      gradff.push_back(invjac * xtensor::xVector<>(-0.5, 0., 0.));
      gradff.push_back(invjac * xtensor::xVector<>(0.5, 0., 0.));
   }

   //  ... --> c est grace a gradloc que gradff est un vecteur de dimension 4 de mVectors
   if (debug)
   {
      for (int j = 0; j < dim + 1; j++)
      {
         printf("node %12.5e %12.5e %12.5e gradient %12.5e %12.5e %12.5e\n", p[j](0), p[j](1), p[j](2), gradff[j](0),
                gradff[j](1), gradff[j](2));
      }
   }
}

// ICI, on suppose qu on connait les valeurs des gradients des FF au point considere
xtensor::xVector<> xElement::getGradInterpoSca(const std::vector<double>& scas)
{
   //  ... calcul du gradient d un champ scalaire interpole avec les fonctions de formes
   std::vector<xtensor::xVector<>> gradloc;
   getGradFF(gradloc);
   return std::inner_product(gradloc.begin(), gradloc.end(), scas.begin(), xtensor::xVector<>(0., 0., 0.));
}

// ICI, on suppose qu on connait les valeurs des gradients des FF au point considere
void xElement::getGradInterpoVec(const std::vector<xtensor::xVector<>>& vecs, xtensor::xTensor2<>& ten)
{
   //  ... calcul du gradient d un champ vectoriel interpole avec les fonctions de formes
   std::vector<xtensor::xVector<>> gradloc;
   getGradFF(gradloc);
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         ten(i, j) = 0.;
         for (int k = 0; k < dim + 1; k++)
         {
            ten(i, j) += vecs[k](i) * gradloc[k](j);
         }
      }
   }

   return;
}

void xElement::xyz2uvw(const xtensor::xPoint& xyz)
{
   xtensor::xVector<> vec(p[0], xyz);
   xtensor::xVector<> res = tinvjac * vec;
   uvw(0) = res(0);
   uvw(1) = res(1);
   uvw(2) = res(2);
   if (type == mEntity::EDGE)
   {
      uvw(0) -= 1.;
   }
   if (type == mEntity::QUAD)
   {
      xmapping::xLagrangeMapping mapping(xmapping::xReferenceElementType::QUAD, p);
      mapping.invert(xyz(0), xyz(1), xyz(2), uvw(0), uvw(1), uvw(2));
   }

   if (type == mEntity::HEX)
   {
      xmapping::xLagrangeMapping mapping(xmapping::xReferenceElementType::HEX, p);
      mapping.invert(xyz(0), xyz(1), xyz(2), uvw(0), uvw(1), uvw(2));
   }
}

//  ... volume de l element
double xElement::getVolume()
{
   if (dim == 0)
   {
      return 1.;
   }
   else if (dim == 1)
   {
      return std::fabs(detjac);
   }
   else if (dim == 2)
   {
      if (type == mEntity::TRI)
      {
         return std::fabs(detjac) / 2;
      }
      else if (type == mEntity::QUAD)
      {
         return std::fabs(detjac);
      }
   }
   else if (dim == 3)
   {
      if (type == mEntity::TET)
      {
         return std::fabs(detjac) / 6;
      }
      else if (type == mEntity::HEX)
      {
         return std::fabs(detjac);
      }
   }
   else
   {
      std::cout << "Warning : getVolume() not coded for the type of element you are trying to use " << std::endl;
      assert(0);
   }

   return 0.;
}

xElementInfo::xElementInfo(mEntity* e) : elem(e)
{
   xElement elem(e);
   volume = elem.getVolume();
   elem.getGradFF(grad_ff_node);
}

xElementInfoManager::xElementInfoPtr xElementInfoManager::getInfo(mEntity* e)
{
   iterator f = infos.find(e);
   if (f != infos.end())
   {
      return f->second;
   }
   xElementInfoPtr info = xElementInfoPtr(new xElementInfo(e));
   infos.insert(value(e, info));
   return info;
}

void xElement::setUvw(const xtensor::xPoint& in) { uvw = in; }
xtensor::xPoint xElement::getUvw() const { return uvw; }

void xElementInfoManager::clear() { infos.clear(); }

}  // namespace xfem
