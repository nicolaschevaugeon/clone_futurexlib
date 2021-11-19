/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include <algorithm>
// mapping
#include "xMapping.h"

namespace xmapping
{
const double edge_vert_pos[2] = {-1., 1.};
const double tri_vert_pos[3][2] = {{0., 0.}, {1., 0.}, {0., 1.}};
const double quad_vert_pos[4][2] = {{-1., -1.}, {1., -1.}, {1., 1.}, {-1., 1.}};
const double tet_vert_pos[4][3] = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
const double prism_vert_pos[8][3] = {{0., 0., -1.}, {1., 0., -1.}, {0., 1., -1.}, {0., 0., 1.}, {1., 0., 1.}, {0., 1., 1.}};
const double hex_vert_pos[8][3] = {{-1., -1., -1.}, {1., -1., -1.}, {1., 1., -1.}, {-1., 1., -1.},
                                   {-1., -1., 1.},  {1., -1., 1.},  {1., 1., 1.},  {-1., 1., 1.}};

const vector3d vertex_normal_of_edge[2] = {{-1., 0., 0.}, {1., 0., 0.}};
const unsigned short edge_vertex_of_tri[3][2] = {{0, 1}, {1, 2}, {2, 0}};
const vector3d edge_normal_of_tri[3] = {{0., -1., 0.}, {1. / sqrt(2.), 1. / sqrt(2.), 0.}, {-1., 0., 0.}};
const unsigned short swapcase_tri[6][3] = {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
const unsigned short edge_vertex_of_quad[4][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}};
const vector3d edge_normal_of_quad[4] = {{0., -1., 0.}, {1., 0., 0.}, {0., 1., 0.}, {-1., 0., 0.}};
const unsigned short swapcase_quad[8][4] = {{0, 1, 2, 3}, {3, 0, 1, 2}, {2, 3, 0, 1}, {1, 2, 3, 0},
                                            {1, 0, 3, 2}, {2, 1, 0, 3}, {3, 2, 1, 0}, {0, 3, 2, 1}};

const unsigned short face_vertex_of_tet[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
const vector3d face_normal_of_tet[4] = {
    {0., 0., -1.}, {0., -1., 0.}, {1. / sqrt(3.), 1. / sqrt(3.), 1. / sqrt(3.)}, {-1., 0., 0.}};

const unsigned short face_vertex_of_hex[6][4] = {{0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
                                                 {3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7}};
const vector3d face_normal_of_hex[6] = {{0., 0., -1.}, {0., -1., 0.}, {1., 0., 0.}, {0., 1., 0.}, {-1., 0., 0.}, {0., 0., 1.}};

std::pair<unsigned short, unsigned short> whichEdge_of_tri(const unsigned short *edge_vertices)
{
   unsigned short ie = 0;
   unsigned short swap = 0;
   for (; ie < 3; ++ie)
      if (std::is_permutation(edge_vertex_of_tri[ie], edge_vertex_of_tri[ie] + 2, edge_vertices)) break;
   if (ie >= 3) throw;
   if (!std::equal(edge_vertex_of_tri[ie], edge_vertex_of_tri[ie] + 2, edge_vertices)) swap = 1;
   return std::make_pair(ie, swap);
}

std::pair<unsigned short, unsigned short> whichEdge_of_quad(const unsigned short *edge_vertices)
{
   unsigned short ie = 0;
   unsigned short swap = 0;
   for (; ie < 4; ++ie)
      if (std::is_permutation(edge_vertex_of_quad[ie], edge_vertex_of_quad[ie] + 2, edge_vertices)) break;
   if (ie >= 4) throw;
   if (!std::equal(edge_vertex_of_quad[ie], edge_vertex_of_quad[ie] + 2, edge_vertices)) swap = 1;
   return std::make_pair(ie, swap);
}

std::pair<unsigned short, unsigned short> whichTri_of_tet(const unsigned short *face_vertices)
{
   unsigned short ie = 0;
   unsigned short swap = 0;
   for (; ie < 4; ++ie)
      if (std::is_permutation(face_vertex_of_tet[ie], face_vertex_of_tet[ie] + 3, face_vertices)) break;
   if (ie >= 4) throw;
   const unsigned short *face_vertex_template = face_vertex_of_tet[ie];
   for (unsigned short swap = 0; swap < 6; ++swap)
   {
      const unsigned short swapped[3] = {face_vertex_template[swapcase_tri[swap][0]], face_vertex_template[swapcase_tri[swap][1]],
                                         face_vertex_template[swapcase_tri[swap][2]]};
      if (std::equal(face_vertex_template, face_vertex_template + 3, swapped)) break;
   }
   if (swap >= 6) throw;
   return std::make_pair(ie, swap);
}

std::pair<unsigned short, unsigned short> whichQuad_of_hex(const unsigned short *face_vertices)
{
   unsigned short ie = 0;
   unsigned short swap = 0;
   for (; ie < 6; ++ie)
      if (std::is_permutation(face_vertex_of_hex[ie], face_vertex_of_hex[ie] + 4, face_vertices)) break;
   if (ie >= 6) throw;
   const unsigned short *face_vertex_template = face_vertex_of_hex[ie];
   for (unsigned short swap = 0; swap < 8; ++swap)
   {
      const unsigned short swapped[4] = {
          face_vertex_template[swapcase_quad[swap][0]],
          face_vertex_template[swapcase_quad[swap][1]],
          face_vertex_template[swapcase_quad[swap][2]],
          face_vertex_template[swapcase_quad[swap][3]],
      };
      if (std::equal(face_vertex_template, face_vertex_template + 4, swapped)) break;
   }
   if (swap >= 8) throw;
   return std::make_pair(ie, swap);
}

void edge_to_tri_map(double s, const unsigned short iedge, const unsigned short swap, double &u, double &v)
{
   if (swap) s = -s;
   const double L0 = (1 - s) / 2;
   const double L1 = 1. - L0;
   unsigned int i0 = edge_vertex_of_tri[iedge][0];
   unsigned int i1 = edge_vertex_of_tri[iedge][1];
   const double *uv0 = tri_vert_pos[i0];
   const double *uv1 = tri_vert_pos[i1];
   u = L0 * uv0[0] + L1 * uv1[0];
   v = L0 * uv0[1] + L1 * uv1[1];
}

void edge_to_quad_map(double s, const unsigned short iedge, const unsigned short swap, double &u, double &v)
{
   if (swap) s = -s;
   const double L0 = (1 - s) / 2;
   const double L1 = 1. - L0;
   const unsigned int i0 = edge_vertex_of_quad[iedge][0];
   const unsigned int i1 = edge_vertex_of_quad[iedge][1];
   const double *uv0 = quad_vert_pos[i0];
   const double *uv1 = quad_vert_pos[i1];
   u = L0 * uv0[0] + L1 * uv1[0];
   v = L0 * uv0[1] + L1 * uv1[1];
}

void face_to_tet_map(const double s, const double t, const unsigned short iface, const unsigned short swap, double &u, double &v,
                     double &w)
{
   const double L_noswap[3] = {1. - s - t, s, t};
   double L[3];
   L[swapcase_tri[swap][0]] = L_noswap[0];
   L[swapcase_tri[swap][1]] = L_noswap[1];
   L[swapcase_tri[swap][2]] = L_noswap[2];
   unsigned int i0 = face_vertex_of_tet[iface][0];
   unsigned int i1 = face_vertex_of_tet[iface][1];
   unsigned int i2 = face_vertex_of_tet[iface][2];
   const double *uvw0 = tet_vert_pos[i0];
   const double *uvw1 = tet_vert_pos[i1];
   const double *uvw2 = tet_vert_pos[i2];
   u = L[0] * uvw0[0] + L[1] * uvw1[0] + L[2] * uvw2[0];
   v = L[0] * uvw0[1] + L[1] * uvw1[1] + L[2] * uvw2[1];
   w = L[0] * uvw0[2] + L[1] * uvw1[2] + L[2] * uvw2[2];
}

void face_to_hex_map(const double s, const double t, const unsigned short iface, const unsigned short swap, double &u, double &v,
                     double &w)
{
   const double L_noswap[4] = {(1. - s) * (1 - t) / 4., (1. + s) * (1 - t) / 4., (1. + s) * (1 + t) / 4.,
                               (1. - s) * (1 + t) / 4.};
   double L[4];
   L[swapcase_quad[swap][0]] = L_noswap[0];
   L[swapcase_quad[swap][1]] = L_noswap[1];
   L[swapcase_quad[swap][2]] = L_noswap[2];
   L[swapcase_quad[swap][3]] = L_noswap[3];

   unsigned int i0 = face_vertex_of_hex[iface][0];
   unsigned int i1 = face_vertex_of_hex[iface][1];
   unsigned int i2 = face_vertex_of_hex[iface][2];
   unsigned int i3 = face_vertex_of_hex[iface][3];

   const double *uvw0 = hex_vert_pos[i0];
   const double *uvw1 = hex_vert_pos[i1];
   const double *uvw2 = hex_vert_pos[i2];
   const double *uvw3 = hex_vert_pos[i3];

   u = L[0] * uvw0[0] + L[1] * uvw1[0] + L[2] * uvw2[0] + L[3] * uvw3[0];
   v = L[0] * uvw0[1] + L[1] * uvw1[1] + L[2] * uvw2[1] + L[3] * uvw3[1];
   w = L[0] * uvw0[2] + L[1] * uvw1[2] + L[2] * uvw2[2] + L[3] * uvw3[2];
}

vector3d compute_normal_to_face_edge(const xMapping &facemapping, const unsigned short *face_vertices, double s)
{
   double u, v;
   vector3d n_uv;
   switch (facemapping.getType())
   {
      case xReferenceElementType::TRI:
      {
         const auto iedge_swap = whichEdge_of_tri(face_vertices);
         const auto iedge = iedge_swap.first;
         const auto swap = iedge_swap.second;
         edge_to_tri_map(s, iedge, swap, u, v);
         n_uv = edge_normal_of_tri[iedge];
         break;
      }
      case xReferenceElementType::QUAD:
      {
         const auto iedge_swap = whichEdge_of_quad(face_vertices);
         const auto iedge = iedge_swap.first;
         const auto swap = iedge_swap.second;
         edge_to_quad_map(s, iedge, swap, u, v);
         n_uv = edge_normal_of_quad[iedge];
         break;
      }
      default:
         throw;
   }
   vector3d n_xyz(n_uv);
   facemapping.pushBack(u, v, 0., 1, &n_xyz);
   return n_xyz.norm();
}

vector3d compute_normal_to_solid_face(const xMapping &solidmapping, const unsigned short *face_vertices, double s, double t)
{
   double u, v, w;
   vector3d n_uvw;
   switch (solidmapping.getType())
   {
      case xReferenceElementType::TET:
      {
         const auto iface_swap = whichTri_of_tet(face_vertices);
         const auto iface = iface_swap.first;
         const auto swap = iface_swap.second;
         face_to_tet_map(s, t, iface, swap, u, v, w);
         n_uvw = face_normal_of_tet[iface];
         break;
      }
      case xReferenceElementType::HEX:
      {
         const auto iface_swap = whichQuad_of_hex(face_vertices);
         const auto iface = iface_swap.first;
         const auto swap = iface_swap.second;
         face_to_hex_map(s, t, iface, swap, u, v, w);
         n_uvw = face_normal_of_hex[iface];
         break;
      }
      default:
         throw;
   }
   vector3d n_xyz(n_uvw);
   solidmapping.pushBack(u, v, w, 1, &n_xyz);
   return n_xyz.norm();
}

xMapping::xMapping(xReferenceElementType _type) : type(_type) {}

xtensor::xPoint xMapping::eval(const xtensor::xPoint &uvw) const
{
   xtensor::xPoint xyz;
   eval(uvw[0], uvw[1], uvw[2], xyz[0], xyz[1], xyz[2]);
   return xyz;
}

tensor xMapping::deval(const xtensor::xPoint &uvw) const
{
   tensor dxdu;
   deval(uvw[0], uvw[1], uvw[2], dxdu(0, 0), dxdu(1, 0), dxdu(2, 0), dxdu(0, 1), dxdu(1, 1), dxdu(2, 1), dxdu(0, 2), dxdu(1, 2),
         dxdu(2, 2));
   return dxdu;
}

xgeom::xBoundingBox xMapping::boundingBox() const
{
   xgeom::xBoundingBox bb;
   boundingBox(bb.min, bb.max);
   return bb;
}

double xMapping::pushBack(double u, double v, double w, size_t vsize, vector3d *gr) const
{
   tensor jInv;
   double detJac = jacInverse(u, v, w, jInv);
   for (size_t i = 0; i < vsize; i++) gr[i] *= jInv;
   return detJac;
}

double xMapping::pushBack(const xtensor::xPoint &uvw, size_t vsize, vector3d *vec) const
{
   return pushBack(uvw[0], uvw[1], uvw[2], vsize, vec);
}

double xMapping::pushBack(double u, double v, double w, size_t vsize, tensor *gr) const
{
   tensor jInv;
   double detJac = jacInverse(u, v, w, jInv);
   for (size_t i = 0; i < vsize; i++) gr[i] = jInv * gr[i];

   /*! \note
       gr[i] is a xTensor2 whose jk element is kth component of vector3d Ni,j
       gr[i][0] = Ni,r,	gr[i][1] = Ni,s,	gr[i][2] = Ni,t
       jInv is the inverse of the Jacobian
           |r,x s,x t,x|
           |r,y s,y t,y|
           |r,x s,z t,z|

       After above matrix multiplication:
       gr[i][0] = Ni,x, gr[i][1] = Ni,y, gr[i][2] = Ni,z
   */

   return detJac;
}

double xMapping::pushBack(const xtensor::xPoint &uvw, size_t vsize, tensor *tens) const
{
   return pushBack(uvw[0], uvw[1], uvw[2], vsize, tens);
}

bool xMapping::invert(double xp, double yp, double zp, double &Upos, double &Vpos, double &Wpos) const
{
#define NR_PRECISION 1.e-6
#define NR_MAX_ITER 50

   double x_est, y_est, z_est;
   double u_new, v_new, w_new;
   double Error = 1.0;
   int iter = 1;
   COG(Upos, Vpos, Wpos);
   tensor InvJacMatrix;
   while (Error > NR_PRECISION && iter < NR_MAX_ITER)
   {
      iter++;
      jacInverse(Upos, Vpos, Wpos, InvJacMatrix);
      eval(Upos, Vpos, Wpos, x_est, y_est, z_est);
      // printf("%f %f %f %f %f %f\n",Upos,Vpos,Wpos,x_est,y_est,z_est);
      u_new = Upos + InvJacMatrix(0, 0) * (xp - x_est) + InvJacMatrix(1, 0) * (yp - y_est) + InvJacMatrix(2, 0) * (zp - z_est);
      v_new = Vpos + InvJacMatrix(0, 1) * (xp - x_est) + InvJacMatrix(1, 1) * (yp - y_est) + InvJacMatrix(2, 1) * (zp - z_est);
      w_new = Wpos + InvJacMatrix(0, 2) * (xp - x_est) + InvJacMatrix(1, 2) * (yp - y_est) + InvJacMatrix(2, 2) * (zp - z_est);
      Error = (u_new - Upos) * (u_new - Upos) + (v_new - Vpos) * (v_new - Vpos) + (w_new - Wpos) * (w_new - Wpos);
      Upos = u_new;
      Vpos = v_new;
      Wpos = w_new;
   }
   if (Error > NR_PRECISION)
   {
      // printf("impossible to find %f %f %f \n",xp,yp,zp);
      return false;
   }
   return true;
}

bool xMapping::invert(const xtensor::xPoint &xyz, xtensor::xPoint &uvw) const
{
   return invert(xyz[0], xyz[1], xyz[2], uvw[0], uvw[1], uvw[2]);
}

// certainly not the fastest way !!
double xMapping::detJac(double u, double v, double w) const
{
   tensor t;
   return jacInverse(u, v, w, t);
}

double xMapping::detJac(const xtensor::xPoint &uvw) const
{
   tensor t;
   return jacInverse(uvw, t);
}

/**
   This function is valid for any mesh mapping
*/
double xMapping::jacInverse(double u, double v, double w, tensor &Invjac) const
{
   double jac[3][3];
   double DetJac = 1;
   deval(u, v, w, jac[0][0], jac[0][1], jac[0][2], jac[1][0], jac[1][1], jac[1][2], jac[2][0], jac[2][1], jac[2][2]);
   switch (type)
   {
      case xReferenceElementType::VERTEX:
         break;
      // may take into account non x-y elements and curved elements, curved lines,...
      case xReferenceElementType::EDGE:
      {
         DetJac = sqrt(jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1] + jac[0][2] * jac[0][2]);
         // vrai is true in french
         int vrai = 0;
         for (int i = 0; i < 3; i++)
         {
            if (jac[0][i] == 0)
            {  // if this component of the normal vector is zero
               vrai = 1;
               for (int j = 0; j < 3; j++)
               {
                  if (j == i)
                     jac[1][j] = 1;  // then the component of the second normal must be one,
                  else
                     jac[1][j] = 0;  // and the other components of the second normal are zero.
               }
               continue;
            }
         }
         // surface equation with normal vector n : n1*X + n2*Y + n3*Z = 0
         if (!vrai)
         {  // looking for the second normal vector in the plane z = 0
            double temp = sqrt(jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1]);
            jac[1][0] = -jac[0][1] / temp;
            jac[1][1] = jac[0][0] / temp;
            jac[1][2] = 0;
         }
         // The third normal vector
         jac[2][0] = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) / DetJac;
         jac[2][1] = (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]) / DetJac;
         jac[2][2] = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) / DetJac;
         break;
      }
      case xReferenceElementType::TRI:
      case xReferenceElementType::QUAD:
      {
         double d3 = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
         double d2 = jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2];
         double d1 = jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1];
         DetJac = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
         jac[2][0] = d1 / DetJac;
         jac[2][1] = d2 / DetJac;
         jac[2][2] = d3 / DetJac;
         break;
      }
      case xReferenceElementType::TET:
      case xReferenceElementType::HEX:
      case xReferenceElementType::PRISM:
      case xReferenceElementType::PYRAMID:
      {
         DetJac = jac[0][0] * jac[1][1] * jac[2][2] + jac[0][2] * jac[1][0] * jac[2][1] + jac[0][1] * jac[1][2] * jac[2][0] -
                  jac[0][2] * jac[1][1] * jac[2][0] - jac[0][0] * jac[1][2] * jac[2][1] - jac[0][1] * jac[1][0] * jac[2][2];
         break;
      }
      default:
      {
         std::cout << "Error :  DetJac not defined for mapping type " << static_cast<int>(type) << " in " << __FILE__ << ":"
                   << __LINE__ << std::endl;
         throw -2345;
      }
   }

   Invjac(0, 0) = (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) / DetJac;
   Invjac(1, 0) = -(jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) / DetJac;
   Invjac(2, 0) = (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) / DetJac;

   Invjac(0, 1) = -(jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1]) / DetJac;
   Invjac(1, 1) = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) / DetJac;
   Invjac(2, 1) = -(jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0]) / DetJac;

   Invjac(0, 2) = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) / DetJac;
   Invjac(1, 2) = -(jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0]) / DetJac;
   Invjac(2, 2) = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) / DetJac;
   return DetJac;
}
double xMapping::jacInverse(const xtensor::xPoint &uvw, tensor &Invjac) const
{
   return jacInverse(uvw[0], uvw[1], uvw[2], Invjac);
}
bool xMapping::inReferenceElement(double u, double v, double w) const
{
   const double eps = 1.e-6;

   switch (type)
   {
      case xReferenceElementType::TRI:
         if (u < -eps || v < -eps || 1. - u - v < -eps) return false;
         break;
      case xReferenceElementType::QUAD:
         if (u < (-1. - eps) || u > (1 + eps) || v < (-1.0 - eps) || v > (1.0 + eps)) return false;
         break;
      case xReferenceElementType::EDGE:
         if (u < -1. - eps || u > 1 + eps) return false;
         break;
      case xReferenceElementType::TET:
         if (u < -eps || v < -eps || w < -eps || 1. - u - v - w < -eps) return false;
         break;
      case xReferenceElementType::HEX:
         if (u < (-1. - eps) || u > (1 + eps) || v < (-1. - eps) || v > (1 + eps) || w < (-1. - eps) || w > (1 + eps))
            return false;
         break;
      case xReferenceElementType::PRISM:
         if (u < -eps || v < -eps || 1. - u - v < -eps || w < (-1. - eps) || w > (1 + eps)) return false;
         break;
      default:
      {
         std::cout << "Error :  inReferenceElement not defined for mapping type " << static_cast<int>(type) << " in " << __FILE__
                   << ":" << __LINE__ << std::endl;
         throw;
      }
   }
   return true;
}

bool xMapping::inReferenceElement(const xtensor::xPoint &uvw) const { return inReferenceElement(uvw[0], uvw[1], uvw[2]); }

bool xMapping::interiorCheck(const xtensor::xPoint &p, xtensor::xPoint &uvw) const
{
   if (!(boundingBox().contains(p))) return false;
   if (!invert(p, uvw)) return false;
   if (!inReferenceElement(uvw)) return false;
   return true;
}

void xMapping::COG(double &u, double &v, double &w) const
{
   switch (type)
   {
      case xReferenceElementType::PRISM:
      case xReferenceElementType::TRI:
         u = v = 0.3333333333333;
         w = 0;
         break;
      case xReferenceElementType::HEX:
      case xReferenceElementType::QUAD:
      case xReferenceElementType::EDGE:
         u = v = w = 0.0;
         break;
      case xReferenceElementType::TET:
         u = v = w = 0.25;
         break;
      default:
      {
         std::cout << "Error : COG unknown for type " << static_cast<int>(type) << " in " << __FILE__ << ":" << __LINE__
                   << std::endl;
         break;  // need to throw an exception here
      }
   }
}

xtensor::xPoint xMapping::COG() const
{
   xtensor::xPoint uvw;
   COG(uvw[0], uvw[1], uvw[2]);
   return uvw;
}

vector3d xMapping::normalVector(const unsigned short *face_vertices, double s, double t) const
{
   switch (type)
   {
      case xReferenceElementType::EDGE:
      {
         unsigned short iv = *face_vertices;
         vector3d n_uvw = vertex_normal_of_edge[iv];
         vector3d n_xyz(n_uvw);
         if (iv)
            pushBack(-1., 0., 0., 1, &n_xyz);
         else
            pushBack(1., 0., 0., 1, &n_xyz);
         return n_xyz.norm();
      }
      case xReferenceElementType::TRI:
      case xReferenceElementType::QUAD:
         return compute_normal_to_face_edge(*this, face_vertices, s);
      case xReferenceElementType::TET:
      case xReferenceElementType::HEX:
         return compute_normal_to_solid_face(*this, face_vertices, s, t);
      default:
         throw;
   }
}
//! compute the exterior normal to the  face face_vertices at point u,v,w  of the face
vector3d xMapping::normalVector(const unsigned short iface, double u, double v, double w) const
{
   vector3d n_uvw;
   switch (type)
   {
      case xReferenceElementType::EDGE:
         n_uvw = vertex_normal_of_edge[iface];
         break;
      case xReferenceElementType::TRI:
         n_uvw = edge_normal_of_tri[iface];
         break;
      case xReferenceElementType::QUAD:
         n_uvw = edge_normal_of_quad[iface];
         break;
      case xReferenceElementType::TET:
         n_uvw = face_normal_of_tet[iface];
         break;
      case xReferenceElementType::HEX:
         n_uvw = face_normal_of_hex[iface];
         break;
      default:
         throw;
   }
   vector3d n_xyz(n_uvw);
   pushBack(u, v, w, 1, &n_xyz);
   return n_xyz.norm();
}

vector3d xMapping::normalVector(double u, double v, double w) const
{
   double jac[3][3];
   deval(u, v, w, jac[0][0], jac[0][1], jac[0][2], jac[1][0], jac[1][1], jac[1][2], jac[2][0], jac[2][1], jac[2][2]);

   switch (type)
   {
      case xReferenceElementType::EDGE:
      {
         const vector3d t(jac[0][0], jac[0][1], jac[0][2]);
         const vector3d z(0, 0, 1);
         vector3d n = t % z;
         return n.norm();
      }
      break;
      case xReferenceElementType::TRI:
      case xReferenceElementType::QUAD:
      {
         vector3d t1(jac[0][0], jac[0][1], jac[0][2]);
         vector3d t2(jac[1][0], jac[1][1], jac[1][2]);
         vector3d n = t1 % t2;
         return n.norm();
      }
      break;
      default:
      {
         std::cout << "Error in xMapping::normalVector : type must be either EDGE, TRI or QUAD" << std::endl;
         throw 1;
      }
   }
}

size_t xMapping::getNbMappingVertices(xReferenceElementType type)
{
   switch (type)
   {
      case xReferenceElementType::VERTEX:
         return 1;
      case xReferenceElementType::EDGE:
         return 2;
      case xReferenceElementType::TRI:
         return 3;
      case xReferenceElementType::QUAD:
         return 4;
      case xReferenceElementType::TET:
         return 4;
      case xReferenceElementType::PYRAMID:
         return 5;
      case xReferenceElementType::PRISM:
         return 6;
      case xReferenceElementType::HEX:
         return 8;
      default:
         throw;
   }
}
xtensor::xPoint xMapping::getMappingVertex(xReferenceElementType type, size_t i)
{
   switch (type)
   {
      case xReferenceElementType::VERTEX:
         assert(i < 1);
         return xtensor::xPoint{0., 0., 0.};
      case xReferenceElementType::EDGE:
         assert(i < 2);
         return xtensor::xPoint{edge_vert_pos[i], 0., 0.};
      case xReferenceElementType::TRI:
         assert(i < 3);
         return xtensor::xPoint{tri_vert_pos[i][0], tri_vert_pos[i][1], 0.};
      case xReferenceElementType::QUAD:
         assert(i < 4);
         return xtensor::xPoint{quad_vert_pos[i][0], quad_vert_pos[i][1], 0.};
      case xReferenceElementType::TET:
         assert(i < 4);
         return xtensor::xPoint{tet_vert_pos[i][0], tet_vert_pos[i][1], tet_vert_pos[i][2]};
      case xReferenceElementType::PRISM:
         assert(i < 6);
         return xtensor::xPoint{prism_vert_pos[i][0], prism_vert_pos[i][1], prism_vert_pos[i][2]};
      case xReferenceElementType::HEX:
         assert(i < 8);
         return xtensor::xPoint{hex_vert_pos[i][0], hex_vert_pos[i][1], hex_vert_pos[i][2]};
      default:
         throw;
   }
}

xtensor::xPoint xMapping::getMappingVertex(size_t i) const { return getMappingVertex(type, i); }

size_t xMapping::getNbMappingVertices() const { return getNbMappingVertices(type); }

}  // namespace xmapping
