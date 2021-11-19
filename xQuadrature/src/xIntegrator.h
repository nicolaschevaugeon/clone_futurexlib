
#ifndef _XQ__INTEGRATOR_H_
#define _XQ__INTEGRATOR_H_
// std
#include <iostream>
// xinterface::aomd
#include "xAOMDEntityUtil.h"
// AOMD
#include "mEntity.h"
// xmapping
#include "xReferenceElement.h"
// xfem
#include "xMapping.h"
// std
#include <iostream>
// xquadrature
#include "xIntPt.h"

namespace xquadrature
{
int getNGQTPts(int);
int getNGQQPts(int);
int getNGQTetPts(int);
int getNGQHPts(int);
xIntPt2d *getGQTPts(int order);
xIntPt2d *getGQQPts(int order);
xIntPt3d *getGQTetPts(int order);
xIntPt3d *getGQHPts(int order);
int GaussLegendre1D(int, double **, double **);

class xIntegrator
{
  public:
   virtual ~xIntegrator() {}
   virtual int nbIntegrationPoints(int order) const = 0;
   virtual void iPoint(int i, int order, double &u, double &v, double &w, double &weight) const = 0;
};

class xGaussIntegrator : public xIntegrator
{
  private:
   xmapping::xReferenceElementType refType;

  public:
   xGaussIntegrator(AOMD::mEntity *e) { refType = xinterface::aomd::mapToxReferenceElementType(e->getType()); };

   xGaussIntegrator(xmapping::xReferenceElementType etype) { refType = etype; };

   virtual int nbIntegrationPoints(int order) const
   {
      switch (refType)
      {
         case xmapping::xReferenceElementType::VERTEX:
            return 1;
         case xmapping::xReferenceElementType::EDGE:
            return order / 2 + 1;  //(order<1)?1:order;
         case xmapping::xReferenceElementType::TRI:
            return (order < 1) ? 1 : getNGQTPts(order);
         case xmapping::xReferenceElementType::TET:
            return (order < 1) ? 1 : getNGQTetPts(order);
         case xmapping::xReferenceElementType::QUAD:
            return getNGQQPts(order);
         case xmapping::xReferenceElementType::HEX:
            return getNGQHPts(order);
         default:
         {
            std::cout << "refType =" << int(refType) << " is not taken into account in xIntegrator L79." << std::endl;
            throw 1;
         }
      }
   };

   virtual void iPoint(int i, int order, double &u, double &v, double &w, double &weight) const
   {
      switch (refType)
      {
         case xmapping::xReferenceElementType::VERTEX:
            u = v = w = 0.0;
            weight = 1.0;
            break;
         case xmapping::xReferenceElementType::EDGE:
         {
            double *pt, *wt;
            GaussLegendre1D(nbIntegrationPoints(order), &pt, &wt);
            u = pt[i];
            v = w = 0.0;
            weight = wt[i];
         }
         break;
         case xmapping::xReferenceElementType::TRI:
         {
            xIntPt2d *pts = getGQTPts(order);
            u = pts[i].pt[0];
            v = pts[i].pt[1];
            w = 0.0;
            weight = 0.5 * pts[i].weight;
         }
         break;
         case xmapping::xReferenceElementType::QUAD:
         {
            xIntPt2d *pts = getGQQPts(order);
            u = pts[i].pt[0];
            v = pts[i].pt[1];
            w = 0.0;
            weight = pts[i].weight;
         }
         break;
         case xmapping::xReferenceElementType::HEX:
         {
            xIntPt3d *pts = getGQHPts(order);
            u = pts[i].pt[0];
            v = pts[i].pt[1];
            w = pts[i].pt[2];
            weight = pts[i].weight;
         }
         break;
         case xmapping::xReferenceElementType::TET:
         {
            const double sixth = 1.0 / 6.0;
            xIntPt3d *pts = getGQTetPts(order);
            u = pts[i].pt[0];
            v = pts[i].pt[1];
            w = pts[i].pt[2];
            weight = sixth * pts[i].weight;
         }
         break;
         default:
         {
            std::cout << "refType =" << int(refType) << " is not taken into account in xIntegrator L 140." << std::endl;
            throw 846;
         }
      }
   };
};

}  // namespace xquadrature
#endif
