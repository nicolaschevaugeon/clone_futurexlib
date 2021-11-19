/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _xMAPPING_
#define _xMAPPING_

#include <vector>
// xtensor
#include "xBoundingBox.h"
#include "xPoint.h"
#include "xTensor2.h"
#include "xVector.h"
// xmapping
#include "xReferenceElement.h"

namespace xmapping
{
using tensor = xtensor::xTensor2<double>;
using vector3d = xtensor::xVector<double>;

/// xMapping is virtual class to define a mapping from a reference domain to the space domain
/*!
 * \brief The xMapping class define a mapping from a reference domain to the space domain.
 * u, v, w refere to coordinate in the reference space, x, y, z refer to coordinate in the physical space.
 * derived class must implement a minima the following pur virtual function :
 *  virtual void eval(double u, double v, double w, double &x, double &y, double &z) const = 0;
 *  virtual void deval(double u, double v, double w, double &dxdu, double &dydu, double &dzdu, double &dxdv, double &dydv,
 *                     double &dzdv, double &dxdw, double &dydw, double &dzdw) const = 0;
 *  virtual void boundingBox(xtensor::xPoint &min, xtensor::xPoint &max) const = 0;
 *  virtual int order() const = 0;
 *  virtual int geomOrder() const = 0;
 */

class xMapping
{
  public:
   xMapping(xReferenceElementType type);
   virtual ~xMapping() = default;
   /// return the type of the refence element
   inline xReferenceElementType getType() const { return type; }
   /// return the coordinates x,y,z in real space corresponding to the local cordinate u,v,w
   virtual void eval(double u, double v, double w, double &x, double &y, double &z) const = 0;
   /// return the coordinates xyz in real space corresponding to the local cordinate uvw
   virtual xtensor::xPoint eval(const xtensor::xPoint &uvw) const;
   /// return the partial derivatives of real coordinates xyz with regard to uvw at local cordinate u,v,w
   virtual void deval(double u, double v, double w, double &dxdu, double &dydu, double &dzdu, double &dxdv, double &dydv,
                      double &dzdv, double &dxdw, double &dydw, double &dzdw) const = 0;
   /// return the partial derivatives of real coordinates xyz with regard to uvw at local cordinate u,v,w  in a tensor such as
   //! (i,j) = dx_i/du_j
   virtual tensor deval(const xtensor::xPoint &uvw) const;

   /// Inversion of the mapping by a newton algorithm. return true if succcessful, false other wise
   virtual bool invert(double x, double y, double z, double &u, double &v, double &w) const;
   /// Inversion of the mapping by a newton algorithm. return true if succcessful, false other wise
   virtual bool invert(const xtensor::xPoint &xyz, xtensor::xPoint &uvw) const;
   /*!
      computation of the inverse of the jacobian at u,v,w (invjac),
      where jacobian = d(x,y,z)/d(u,v,w) is a 3x3 matrix
      returns the determinant of the jacobian (not of the inverse).
    */
   virtual double jacInverse(double u, double v, double w, tensor &invjac) const;
   /*!
      computation of the inverse of the jacobian at point uvw (invjac),
      where jacobian = d(x,y,z)/d(u,v,w) is a 3x3 matrix
      returns the determinant of the jacobian (not of the inverse).
    */
   virtual double jacInverse(const xtensor::xPoint &uvw, tensor &invjac) const;
   //!    returns determinant of jacobian at (u,v,w)
   virtual double detJac(double u, double v, double w) const;
   //!    returns determinant of jacobian at point uvw
   virtual double detJac(const xtensor::xPoint &uvw) const;

   /**
      tells if a point (u,v,w) is inside the reference element or not
  For EDGE, return true if ( -1.-eps <= u <= 1.+eps(
  for TRI, return true if (u >= - && v >= -eps && 1.-u-v => -eps)
  for QUAD, return true if (u >= (-1.-eps) && u <= (1+eps) && v >= (-1.0-eps) && v <= (1.0+eps))
  for TET, return true if  (u => (-eps && v => -eps && w => -eps && 1.-u-v-w => -eps)
  for HEX, return true if  ( (1+eps) => u, v, w => (-1.-eps)==)
  for PRISM, return true if (-1-eps <= w <= 1+eps && u,v >= -eps && 1.-u-v >= -eps )
   */
   virtual bool inReferenceElement(double u, double v, double w) const;
   virtual bool inReferenceElement(const xtensor::xPoint &uvw) const;
   /**
      checks if a point in REAL coordinates is inside the element ent
      if it's the case, returns local coordinates
    */
   virtual bool interiorCheck(const xtensor::xPoint &p, xtensor::xPoint &uvw) const;
   ///   returns the local coordinates of the COG of the reference element
   virtual void COG(double &u, double &v, double &w) const;
   ///   returns the local coordinates of the COG of the reference element
   virtual xtensor::xPoint COG() const;
   /// Computes gradients at local cordinate u,v,w real world knowing them in reference world, vector version
   //! vsize is the number of vector to push,
   //! on entry vec is a pointer to contiguous vector3d, supposed to represent partialderivative with regrad to u,v,w.
   //! on exit  vec is a pointer to contiguous vector3d, supposed to represent gradient in x,y,z.
   //! the returned value is determinent of the jacobian
   //! \note could use span when switching to c++20
   virtual double pushBack(double u, double v, double w, size_t vsize, vector3d *vec) const;
   /// Computes gradients at local cordinate point uvw real world knowing them in reference world, vector version
   //! vsize is the number of vector to push,
   //! on entry vec is a pointer to contiguous vector3d, supposed to represent partialderivative with regrad to u,v,w.
   //! on exit  vec is a pointer to contiguous vector3d, supposed to represent gradient in x,y,z.
   //! the returned value is determinent of the jacobian
   virtual double pushBack(const xtensor::xPoint &uvw, size_t vsize, vector3d *vec) const;
   inline double pushBack(double u, double v, double w, std::vector<vector3d> &vec) const
   {
      return pushBack(u, v, w, vec.size(), vec.data());
   }
   inline double pushBack(const xtensor::xPoint &uvw, std::vector<vector3d> &vec) const
   {
      return pushBack(uvw, vec.size(), vec.data());
   }
   /// Computes gradients in real world knowing them in reference world, tensor version
   virtual double pushBack(double u, double v, double w, size_t vsize, tensor *tens) const;
   /// Computes gradients in real world knowing them in reference world, tensor version
   virtual double pushBack(const xtensor::xPoint &uvw, size_t vsize, tensor *tens) const;

   /// compute the exterior normal to the  face face_vertices at point s,t of the face
   /*!
    *  compute the normal to the face (resp line) given by face_vertices
    * in term of the nodes number of the reference solid (resp face) upon with the mapping is constructed.
    *  s,t are the local parameter in the reference space of face (resp s in the reference space of the line).
    * It proceed by finding out if the face is really a valid face of the solid,
    * if so which one and how it is rotated compared to the corresponding reference face on the solid.
    * Then it compute the coordinate u,v,w in the solid reference face corresponding to point s,t on the face, taking into account
    * face number and rotations. It then select the normal in refence space according to the face number, and push it back to
    * geometric space.
    */
   virtual vector3d normalVector(const unsigned short *face_vertices, double s, double t) const;

   /// compute the exterior normal to the  numbered iface face relative to the underlying solid at point u,v,w  of the solid.
   /*!
    *  iface is the relative face number to the solid we want to take the normal of. From this value, the normal in ref space is
    * retrived. it is then push to geometric space using the jacobian at u,v,w. It also work on face mapping, then iface must
    * refer to the edge number. It also work on line mapping : then iface must refer to the vertex mapping.
    */
   virtual vector3d normalVector(const unsigned short iface, double u, double v, double w) const;

   ///  computes the normal vector to a face or edge.
   /*!
    * compute the "natural" normal to a face mapping by computing dxdu ^ dxdv.
    *  Direction of the normal of course depends on the order of the nodes in the face up to a sign.
    * it also work with line, making the assumption that the normal is orthogonal to the third axis (dxdv set to e3).
    */
   virtual vector3d normalVector(double u, double v, double w) const;
   /// computes the bounding box of the mapping in real space.
   virtual void boundingBox(xtensor::xPoint &min, xtensor::xPoint &max) const = 0;
   /// computes the bounding box of the mapping in real space.
   virtual xgeom::xBoundingBox boundingBox() const;
   /// return the polynomial order of the mapping (if it make sense ...), used for chosing integration rule
   virtual int order() const = 0;
   /// return the polynomial order of the mapping (if it make sense ...), used for chosing integration rule
   //! similar to order. probably here for historical reason ... \note I need to figure this out
   virtual int geomOrder() const = 0;
   size_t getNbMappingVertices() const;
   xtensor::xPoint getMappingVertex(size_t i) const;
   static size_t getNbMappingVertices(xReferenceElementType type);
   static xtensor::xPoint getMappingVertex(xReferenceElementType type, size_t i);

  protected:
   xReferenceElementType type;
};

}  // namespace xmapping
#endif
