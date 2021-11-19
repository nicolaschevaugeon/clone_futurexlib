// xcut
#include "xPhysSurfVLS.h"

#include "xPhysSurfByTagging.h"
// xtensor
#include "xPoint.h"
// xfem
#include "xMesh.h"
// AOMD
#include "mAOMD.h"
#include "mFace.h"
#include "mVertex.h"
// boost
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace xcut
{
struct delete_ptr
{
   template <typename T>
   void operator()(T *pPtr)
   {
      delete pPtr;
   }
};
// tolerance to accept cut in invalid order as valid
const double eps2 = 1.e-7;

// tolerance to relative volume of a element as non nul to be considered as valid and in the support
const double eps4 = 1.e-2;

bool edgeCut(const xtensor::xPoint &p0, const xtensor::xPoint &p1, const xtensor::xVector<> &v, int which, double &s)
{
   const double eps = 1.e-7;
   const double n = v.mag();
   xtensor::xVector<> p0p1(p0, p1);
   const double m = p0p1.mag();
   if (n < eps)
   {
      s = (which == 0) ? 0. : 1.;
      return true;
   }
   if (m < eps) throw -1;
   const xtensor::xPoint pOnPlane = (which == 0) ? p0 + v : p1 + v;
   const double a = v(0) / n;
   const double b = v(1) / n;
   const double c = v(2) / n;
   const double d = a * pOnPlane(0) + b * pOnPlane(1) + c * pOnPlane(2);
   // Linear equation for s0, parameter of the intersection point
   // es = f.
   const double e = (p0p1(0) * a + p0p1(1) * b + p0p1(2) * c) / m;
   if (fabs(e) > eps)
   {
      const double f = (d - p0(0) * a - p0(1) * b - p0(2) * c) / m;
      s = f / e;
      return true;
   }
   return false;
}

void edgeCut(const xtensor::xPoint &p0, const xtensor::xPoint &p1, const xfem::xVectorLevelSetData *v0data,
             const xfem::xVectorLevelSetData *v1data, std::vector<xtensor::xPoint> &cutpoints)
{
   const double eps = 1.e-6;  // 0.; //1.e-2;

   cutpoints.clear();
   double s0 = -2.;
   double s1 = -2.;
   bool cf0;
   bool cf1;
   if (v0data && v1data && v0data->doExist() && v1data->doExist())
   {
      int inout0 = v0data->getInOut();
      int inout1 = v1data->getInOut();
      xtensor::xVector<> v0 = v0data->getDirectionToIso();
      xtensor::xVector<> v1 = v1data->getDirectionToIso();
      if ((inout0 + inout1) < 0) return;  // case -- -0 0- 00 : no cut.
      if ((inout0 == 1) && (inout1 != 1))
      {  // case +0 or +- : 1 cut, priority from +, - or 0, last resort ls cut.
         cf0 = edgeCut(p0, p1, v0, 0, s0);
         if (cf0 && s0 >= 0 && s0 <= 1)
         {
            if (s0 <= eps) s0 = eps;
            if (s0 >= (1 - eps)) s0 = 1. - eps;
            cutpoints.push_back(p0 * (1 - s0) + p1 * (s0));
            return;
         }
         else if (inout1 == 0)
         {
            double s = 1. - eps;
            //   double s =1.;
            cutpoints.push_back(p0 * (1 - s) + p1 * (s));
            // cutpoints.push_back(p1);
            return;
         }
         else
         {  // case +-
            cf1 = edgeCut(p0, p1, v1, 1, s1);
            if (cf1 && s1 >= 0 && s1 <= 1)
            {  // a cut is found
               if (s1 <= eps) s1 = eps;
               if (s1 >= (1 - eps)) s1 = 1. - eps;
               cutpoints.push_back(p0 * (1 - s1) + p1 * (s1));
               return;
            }
            else
            {  // no cut found from 1 : last resort, cut using vls.
               // double s = 1./(-v1.mag() - v0.mag());
               double s = v0.mag() / (v1.mag() + v0.mag());
               if (s <= eps) s = eps;
               if (s >= (1 - eps)) s = 1. - eps;
               cutpoints.push_back(p0 * (1. - s) + p1 * (s));
               return;
            }
         }
      }
      if ((inout0 != 1) && (inout1 == 1))
      {  // case 0+ or -+ : 1 cut, priority from +, - or 0, last resort ls cut.
         cf1 = edgeCut(p0, p1, v1, 1, s1);
         if (cf1 && s1 >= 0 && s1 <= 1)
         {
            if (s1 <= eps) s1 = eps;
            if (s1 >= (1 - eps)) s1 = 1. - eps;
            cutpoints.push_back(p0 * (1. - s1) + p1 * (s1));
            return;
         }
         else if (inout0 == 0)
         {
            double s = eps;
            cutpoints.push_back(p0 * (1. - s) + p1 * (s));
            return;
         }
         else
         {  // case -+
            cf0 = edgeCut(p0, p1, v0, 0, s0);
            if (cf0 && s0 >= 0 && s0 <= 1)
            {  // a cut is found
               if (s0 <= eps) s0 = eps;
               if (s0 >= (1 - eps)) s0 = 1. - eps;
               cutpoints.push_back(p0 * (1. - s0) + p1 * (s0));
               return;
            }
            else
            {  // no cut found from 1 : last resort, cut using vls.
               double s = v0.mag() / (v1.mag() + v0.mag());
               if (s <= eps) s = eps;
               if (s >= (1 - eps)) s = 1. - eps;
               cutpoints.push_back(p0 * (1. - s) + p1 * (s));
               return;
            }
         }
      }
      cf0 = edgeCut(p0, p1, v0, 0, s0);
      cf1 = edgeCut(p0, p1, v1, 1, s1);
      bool c0v = (cf0 && s0 >= 0 && s0 <= 1);
      bool c1v = (cf1 && s1 >= 0 && s1 <= 1);
      if (c0v && c1v)
      {  // 2 cut have been found.
         if (s0 > s1)
         {  // cut in invalid order :  the two cut are at the same point
            if (s0 - s1 > eps2) return;
            s0 = (s0 + s1) / 2.;
            s1 = s0;
         }
         cutpoints.push_back(p0 * (1. - s0) + p1 * (s0));
         cutpoints.push_back(p0 * (1. - s1) + p1 * (s1));
         return;
      }
      else if (c0v)
      {  // only cut from 1 is valide : two cut at position s0
         //   cutpoints.push_back(p0*(1.-s0) + p1*(s0) ) ;
         // cutpoints.push_back(p0*(1.-s0) + p1*(s0) ) ;
         return;
      }
      else if (c1v)
      {
         // cutpoints.push_back(p0*(1.-s1) + p1*(s1) ) ;
         // cutpoints.push_back(p0*(1.-s1) + p1*(s1) ) ;
         return;
      }
      return;
   }
   // }
   // if one point at least is close to the skeleton
   else if (v0data && v1data)
   {
      // if both are close to the skeleton no cut
      if (!v0data->doExist() && !v1data->doExist()) return;
      // if v0 is close to the skeleton but v1 not
      // => inout0=1
      if (!v0data->doExist() && v1data->doExist())
      {
         // int inout0 = v0data->getInOut();
         assert(v0data->getInOut() == 1);
         const int inout1 = v1data->getInOut();
         // case +s+
         if (inout1 == 1)
         {
            // no cut
            return;
         }
         // case +s0
         else if (inout1 == 0)
         {
            cutpoints.push_back(p1);
            return;
         }
         // case +s-
         else
         {
            xtensor::xVector<> v1 = v1data->getDirectionToIso();
            cf1 = edgeCut(p0, p1, v1, 1, s1);
            if (cf1 && s1 >= 0 && s1 <= 1)  // a cut is found
            {
               if (s1 <= eps) s1 = eps;
               if (s1 >= (1 - eps)) s1 = 1. - eps;
               cutpoints.push_back(p0 * (1 - s1) + p1 * (s1));
               return;
            }
            else
            {  // no cut found from 1 : last resort, cut using  vls radius.
               xtensor::xVector<> edge_vect(p0, p1);
               assert(edge_vect.mag() > 0.);
               double s = v1.mag() / edge_vect.mag();
               if (s <= eps) s = eps;
               if (s >= (1. - eps)) s = 1. - eps;
               if (s > 1.)
               {
                  cout << "Warning : case +s- with no vls cut : last resort use a radius  greater then edge lenght; s=" << s
                       << endl;
               }
               cutpoints.push_back(p0 * s + p1 * (1. - s));
               return;
            }
         }
      }
      // if v1 is close to the skeleton but v0 not
      // => inout1=1
      else
      {
         int inout0 = v0data->getInOut();

         assert(v1data->getInOut());
         // case +s+
         if (inout0 == 1)
         {
            // no cut
            return;
         }
         // case 0+s
         else if (inout0 == 0)
         {
            cutpoints.push_back(p0);
            return;
         }
         // case -+s
         else
         {
            xtensor::xVector<> v0 = v0data->getDirectionToIso();
            cf0 = edgeCut(p0, p1, v0, 0, s0);
            if (cf0 && s0 >= 0 && s0 <= 1)
            {  // a cut is found
               if (s0 <= eps) s0 = eps;
               if (s0 >= (1 - eps)) s0 = 1. - eps;
               cutpoints.push_back(p0 * (1. - s0) + p1 * (s0));
               return;
            }
            else
            {  // no cut found from 0 : last resort, cut using  vls radius.
               xtensor::xVector<> edge_vect(p0, p1);
               double s = v0.mag() / edge_vect.mag();
               if (s <= eps) s = eps;
               if (s >= (1. - eps)) s = 1. - eps;
               if (s > 1.)
               {
                  cout << "Warning : case +s- with no vls cut : last resort use a radius  greater then edge lenght; s=" << s
                       << endl;
               }
               cutpoints.push_back(p0 * (1. - s) + p1 * (s));
               return;
            }
         }
      }
   }

   if (!(v0data && v1data))
   {
      cout << "WARNING : no data on both v0 and V1" << v0data << " " << v1data << endl;
      return;
   }
   //    cout << "arg !" << v0data << " " << v1data  << endl;
}

AOMD::mEntity *new_mEntityFromTriangle(const triangle &in, pGEntity classification)
{
   return new AOMD::mFace(in[0], in[1], in[2], classification);
}

void cutFacePhysSurfVLS(AOMD::mEntity *face, xfem::xMesh &boundarymesh,
                        const std::vector<const xfem::xVectorLevelSetData *> &nodesVLSdata,
                        const std::array<xVLSEdgeCutData *, 3> &edgesVLSdata, xVLSTriangleCutData &tcutdata)
{
   const double eps2 = 1.e-6;
   tcutdata.clear();
   AOMD::mMesh &boundarymmesh = boundarymesh.getMesh();
   pGEntity GEn1 = boundarymmesh.getGEntity(1, 1);
   // pGEntity Gface  =  face->getClassification();
   size_t nbnodes = size_t(face->size(0));
   if (nbnodes != 3) throw;
   int nbout = 0;
   int nbskel = 0;
   for (size_t i = 0; i < nbnodes; ++i)
   {
      if ((nodesVLSdata[i]->getInOut() == 1)) nbout++;
      if (!nodesVLSdata[i]->doExist()) nbskel++;
   }

   if (nbskel) switch (nbskel)
      {
         case 1:
         {
            if (nbout == 1)
            {
               /* case :                     front :
               //           +s
               //          / \
               //         /   \
               //        *-----*             *-----*
               //       /       \
               //    -/0 ________ -/0                 */
               size_t iout = 0;
               while (nodesVLSdata[iout]->getInOut() != 1 || !nodesVLSdata[iout]->doExist()) ++iout;
               const size_t ie0 = (3 + iout) % 3;
               const size_t ie1 = (2 + iout) % 3;
               xVLSEdgeCutData *vlsedgecutdata0 = edgesVLSdata[ie0];
               xVLSEdgeCutData *vlsedgecutdata1 = edgesVLSdata[ie1];
               assert((vlsedgecutdata0->nbcuts() == 1) && (vlsedgecutdata1->nbcuts() == 1));
               AOMD::mVertex *v0 = (*vlsedgecutdata0)(0);
               AOMD::mVertex *v1 = (*vlsedgecutdata1)(0);

               boundarymmesh.createEdge(v0, v1, GEn1);
               AOMD::mVertex *vTout = static_cast<AOMD::mVertex *>(face->get(0, int(iout)));
               AOMD::mVertex *vTin0 = static_cast<AOMD::mVertex *>(face->get(0, int((iout + 1) % 3)));
               AOMD::mVertex *vTin1 = static_cast<AOMD::mVertex *>(face->get(0, int((iout + 2) % 3)));
               triangle TriangleOut = {vTout, v0, v1};
               triangle TriangleIn0 = {v0, v1, vTin0};
               triangle TriangleIn1 = {v1, vTin0, vTin1};
               tcutdata.TrianglesOut = {TriangleOut};
               tcutdata.TrianglesIn = {TriangleIn0, TriangleIn1};
               break;
            }
            else if (nbout == 3)
            {
               // 2+ 1+s : no cut
               tcutdata.fromMeshOut = face;
               break;
            }
            // else nbout==2 => same as 2+s and nbout=2 => pass to it
         }
         case 2:
         {
            if (nbout == 2)
            {
               int iedge = -1;
               for (int i = 0; i < 3; ++i)
               {
                  AOMD::mEntity *edge = face->get(1, i);
                  AOMD::mVertex *v0 = static_cast<AOMD::mVertex *>(edge->get(0, 0));
                  AOMD::mVertex *v1 = static_cast<AOMD::mVertex *>(edge->get(0, 1));
                  const size_t iv0 = (face->get(0, 0) == v0) ? 0 : ((face->get(0, 1) == v0) ? 1 : 2);
                  const size_t iv1 = (face->get(0, 0) == v1) ? 0 : ((face->get(0, 1) == v1) ? 1 : 2);
                  const xfem::xVectorLevelSetData *nodesVLSdata0 = nodesVLSdata[iv0];
                  const xfem::xVectorLevelSetData *nodesVLSdata1 = nodesVLSdata[iv1];
                  if ((nodesVLSdata0->getInOut() == 1) && (nodesVLSdata1->getInOut() == 1))
                  {
                     iedge = i;
                     break;
                  }
               }
               if (iedge == -1) throw;
               xVLSEdgeCutData *vlsedgecutdatainin = edgesVLSdata[size_t(iedge)];
               if (!vlsedgecutdatainin)
               {
                  /* case :
                  //          -/0           front :
                  //          / \
                  //         /   \
                  //        *-----*         *-----*
                  //       /       \
                  //      +s_______ +s             */
                  xVLSEdgeCutData *vlsedgecutdata0 = edgesVLSdata[(iedge + 2) % 3];
                  xVLSEdgeCutData *vlsedgecutdata1 = edgesVLSdata[(iedge + 1) % 3];
                  assert(vlsedgecutdata0 && vlsedgecutdata1);
                  assert(vlsedgecutdata0->nbcuts() == 1 && vlsedgecutdata1->nbcuts() == 1);
                  AOMD::mVertex *v0 = (*vlsedgecutdata0)(0);
                  AOMD::mVertex *v1 = (*vlsedgecutdata1)(0);
                  boundarymmesh.createEdge(v0, v1, GEn1);
                  AOMD::mVertex *vTin = static_cast<AOMD::mVertex *>(face->get(0, (iedge + 2) % 3));
                  AOMD::mVertex *vTout0 = static_cast<AOMD::mVertex *>(face->get(0, iedge));
                  AOMD::mVertex *vTout1 = static_cast<AOMD::mVertex *>(face->get(0, (iedge + 1) % 3));
                  triangle TriangleIn = {vTin, v0, v1};
                  triangle TriangleOut0 = {v0, v1, vTout0};
                  triangle TriangleOut1 = {v1, vTout0, vTout1};
                  tcutdata.TrianglesIn = {TriangleIn};
                  tcutdata.TrianglesOut = {TriangleOut0, TriangleOut1};
               }
               else
               {
                  throw -1;
               }
            }
            else
            {
               tcutdata.fromMeshOut = face;
            }
            break;
         }
         case 3:
         {
            tcutdata.fromMeshOut = face;
            break;
         }

         default:
            throw;
      }
   else
      switch (nbout)
      {
         case 0:
         {  // nocut
            tcutdata.fromMeshIn = face;
            break;
         }
         case 1:
         {
            /* case :                     front :
            //           +
            //          / \
            //         /   \
            //        *-----*             *-----*
            //       /       \
            //    -/0 ________ -/0                */
            size_t iout = 0;
            while (nodesVLSdata[iout]->getInOut() != 1) ++iout;
            xVLSEdgeCutData *vlsedgecutdata0 = edgesVLSdata[(3 + iout) % 3];
            xVLSEdgeCutData *vlsedgecutdata1 = edgesVLSdata[(2 + iout) % 3];
            assert(vlsedgecutdata0 && vlsedgecutdata1);
            assert((vlsedgecutdata0->nbcuts() == 1) && (vlsedgecutdata1->nbcuts() == 1));
            AOMD::mVertex *v0 = (*vlsedgecutdata0)(0);
            AOMD::mVertex *v1 = (*vlsedgecutdata1)(0);
            boundarymmesh.createEdge(v0, v1, GEn1);
            AOMD::mVertex *vTout = static_cast<AOMD::mVertex *>(face->get(0, int(iout)));
            AOMD::mVertex *vTin0 = static_cast<AOMD::mVertex *>(face->get(0, int((iout + 1) % 3)));
            AOMD::mVertex *vTin1 = static_cast<AOMD::mVertex *>(face->get(0, int((iout + 2) % 3)));
            triangle TriangleOut = {vTout, v0, v1};
            triangle TriangleIn0 = {v0, v1, vTin0};
            triangle TriangleIn1 = {v1, vTin0, vTin1};
            tcutdata.TrianglesOut = {TriangleOut};
            tcutdata.TrianglesIn = {TriangleIn0, TriangleIn1};
            break;
         }
         case 2:
         {
            int iedge = -1;
            for (int i = 0; i < 3; ++i)
            {
               AOMD::mEntity *edge = face->get(1, i);
               AOMD::mVertex *v0 = static_cast<AOMD::mVertex *>(edge->get(0, 0));
               AOMD::mVertex *v1 = static_cast<AOMD::mVertex *>(edge->get(0, 1));
               const size_t iv0 = (face->get(0, 0) == v0) ? 0 : ((face->get(0, 1) == v0) ? 1 : 2);
               const size_t iv1 = (face->get(0, 0) == v1) ? 0 : ((face->get(0, 1) == v1) ? 1 : 2);
               const xfem::xVectorLevelSetData *nodesVLSdata0 = nodesVLSdata[iv0];
               const xfem::xVectorLevelSetData *nodesVLSdata1 = nodesVLSdata[iv1];

               if ((nodesVLSdata0->getInOut() == 1) && (nodesVLSdata1->getInOut() == 1))
               {
                  iedge = i;
                  break;
               }
            }
            if (iedge == -1) throw;
            xVLSEdgeCutData *vlsedgecutdatainin = edgesVLSdata[size_t(iedge)];
            if (!vlsedgecutdatainin)
            {
               /* case :
               //          -/0           front :
               //          / \
               //         /   \
               //        *-----*         *-----*
               //       /       \
               //      + _______ +                 */
               xVLSEdgeCutData *vlsedgecutdata0 = edgesVLSdata[(iedge + 2) % 3];
               xVLSEdgeCutData *vlsedgecutdata1 = edgesVLSdata[(iedge + 1) % 3];
               assert(vlsedgecutdata0 && vlsedgecutdata1);
               assert(vlsedgecutdata0->nbcuts() == 1 && vlsedgecutdata1->nbcuts() == 1);
               AOMD::mVertex *v0 = (*vlsedgecutdata0)(0);
               AOMD::mVertex *v1 = (*vlsedgecutdata1)(0);
               boundarymmesh.createEdge(v0, v1, GEn1);
               AOMD::mVertex *vTin = static_cast<AOMD::mVertex *>(face->get(0, (iedge + 2) % 3));
               AOMD::mVertex *vTout0 = static_cast<AOMD::mVertex *>(face->get(0, iedge));
               AOMD::mVertex *vTout1 = static_cast<AOMD::mVertex *>(face->get(0, (iedge + 1) % 3));
               triangle TriangleIn = {vTin, v0, v1};
               triangle TriangleOut0 = {v0, v1, vTout0};
               triangle TriangleOut1 = {v1, vTout0, vTout1};
               tcutdata.TrianglesIn = {TriangleIn};
               tcutdata.TrianglesOut = {TriangleOut0, TriangleOut1};
            }
            else
            {
               /* case :
               //          -/0
               //          / \            front :
               //         /   \
               //        *     *           *      *
               //       / \   / \           \    /
               //      + -*--*-- +           *  *
               */

               size_t iedgeL = (iedge + 2) % 3;
               size_t iedgeR = (iedge + 1) % 3;
               AOMD::mVertex *vToutL = static_cast<AOMD::mVertex *>(face->get(0, iedge));
               AOMD::mVertex *vToutR = static_cast<AOMD::mVertex *>(face->get(0, (iedge + 1) % 3));
               AOMD::mVertex *vTin = static_cast<AOMD::mVertex *>(face->get(0, (iedge + 2) % 3));
               if (face->get(1, iedge)->get(0, 0) != face->get(0, iedge))
               {
                  swap(iedgeL, iedgeR);
                  swap(vToutL, vToutR);
               }
               assert(vlsedgecutdatainin->nbcuts() == 2);
               AOMD::mVertex *voniniL = (*vlsedgecutdatainin)(0);
               AOMD::mVertex *voniniR = (*vlsedgecutdatainin)(1);
               xVLSEdgeCutData *vlsedgecutdataL = edgesVLSdata[iedgeL];
               xVLSEdgeCutData *vlsedgecutdataR = edgesVLSdata[iedgeR];
               assert((vlsedgecutdataL->nbcuts() == 1) && (vlsedgecutdataR->nbcuts() == 1));
               AOMD::mVertex *vL = (*vlsedgecutdataL)(0);
               AOMD::mVertex *vR = (*vlsedgecutdataR)(0);
               boundarymmesh.createEdge(voniniL, vL, GEn1);
               boundarymmesh.createEdge(voniniR, vR, GEn1);
               triangle TriangleIn0 = {vTin, vL, vR};
               triangle TriangleIn1 = {vL, voniniL, vR};
               triangle TriangleIn2 = {vR, voniniL, voniniR};
               triangle TriangleOut0 = {vL, vToutL, voniniL};
               triangle TriangleOut1 = {vR, voniniR, vToutR};
               tcutdata.TrianglesIn = {TriangleIn0, TriangleIn1, TriangleIn2};
               tcutdata.TrianglesOut = {TriangleOut0, TriangleOut1};
            }
            break;
         }
         case 3:
         {
            int nbedgecut = 0;
            for (size_t i = 0; i < 3; ++i)
               if (edgesVLSdata[i]) ++nbedgecut;
            switch (nbedgecut)
            {
               case 0:
               {
                  tcutdata.fromMeshOut = face;
                  break;
               }
               case 1:
               {
                  /* case :
                  //           +A                 front :
                  //          / \
                  //       v0*   \	         *
                  //        / \   \	        /
                  //    v1 * \  \  \         *
                  //      /     \ \ \
                  //     + --------  +C
                  //     B
                  */
                  size_t iedgecut = 0;
                  while (edgesVLSdata[iedgecut] == nullptr) ++iedgecut;
                  const xVLSEdgeCutData &edgecutdata = *edgesVLSdata[iedgecut];
                  assert(edgecutdata.nbcuts() == 2);
                  AOMD::mVertex *v0 = edgecutdata(0);
                  AOMD::mVertex *v1 = edgecutdata(1);
                  AOMD::mVertex *A = static_cast<AOMD::mVertex *>(face->get(1, int(iedgecut))->get(0, 0));
                  AOMD::mVertex *B = static_cast<AOMD::mVertex *>(face->get(1, int(iedgecut))->get(0, 1));
                  AOMD::mVertex *C = static_cast<AOMD::mVertex *>(face->get(0, int((iedgecut + 2) % 3)));
                  xtensor::xPoint p0 = v0->point();
                  xtensor::xPoint p1 = v1->point();
                  xtensor::xVector<> p0p1(p0, p1);
                  const xtensor::xPoint pA = A->point();
                  const xtensor::xPoint pB = B->point();
                  xtensor::xVector<> pApB(pA, pB);
                  triangle TriangleOut0 = {v0, C, A};

                  if (p0p1.mag() / pApB.mag() <= eps2)
                  {
                     triangle TriangleOut1 = {v0, B, C};
                     tcutdata.TrianglesOut = {TriangleOut0, TriangleOut1};
                  }
                  else
                  {
                     triangle TriangleOut1 = {v0, v1, C};
                     triangle TriangleOut2 = {v1, B, C};
                     tcutdata.TrianglesOut = {TriangleOut0, TriangleOut1, TriangleOut2};
                  }
                  boundarymmesh.createEdge(v0, v1, GEn1);
                  break;
               }
               case 2:
               {
                  /* case :
                  //           +         front :
                  //          / \
                  //         *-- *         *-- *
                  //        /     \
                  //       *------ *     *------ *
                  //      /         \
                  //     + --------  +
                  */
                  size_t iedgenotcut = 0;
                  while (edgesVLSdata[iedgenotcut] != nullptr) ++iedgenotcut;
                  AOMD::mVertex *vTdown0 = static_cast<AOMD::mVertex *>(face->get(0, int(iedgenotcut)));
                  AOMD::mVertex *vTdown1 = static_cast<AOMD::mVertex *>(face->get(0, int((iedgenotcut + 1) % 3)));
                  AOMD::mVertex *vTup = static_cast<AOMD::mVertex *>(face->get(0, int((iedgenotcut + 2) % 3)));
                  const size_t iedge0 = (iedgenotcut + 2) % 3;
                  const size_t iedge1 = (iedgenotcut + 1) % 3;
                  AOMD::mEntity *edge0 = face->get(1, int(iedge0));
                  AOMD::mEntity *edge1 = face->get(1, int(iedge1));
                  xVLSEdgeCutData *edgecutdata0 = edgesVLSdata[iedge0];
                  xVLSEdgeCutData *edgecutdata1 = edgesVLSdata[iedge1];
                  assert((edgecutdata0->nbcuts() == 2) && (edgecutdata1->nbcuts() == 2));
                  AOMD::mVertex *v0up = (*edgecutdata0)(0);
                  AOMD::mVertex *v0down = (*edgecutdata0)(1);
                  AOMD::mVertex *v1up = (*edgecutdata1)(0);
                  AOMD::mVertex *v1down = (*edgecutdata1)(1);
                  if (edge0->get(0, 0) != vTup) swap(v0up, v0down);
                  if (edge1->get(0, 0) != vTup) swap(v1up, v1down);
                  boundarymmesh.createEdge(v0up, v1up, GEn1);
                  boundarymmesh.createEdge(v0down, v1down, GEn1);
                  triangle TOup = {vTup, v0up, v1up};
                  triangle Tin0 = {v0up, v0down, v1up};
                  triangle Tin1 = {v1up, v0down, v1down};
                  triangle TOdown0 = {v0down, vTdown0, v1down};
                  triangle TOdown1 = {v1down, vTdown0, vTdown1};
                  tcutdata.TrianglesOut = {TOup, TOdown0, TOdown1};
                  tcutdata.TrianglesIn = {Tin0, Tin1};
                  break;
               }
               case 3:
               {
                  /* case :
                  //           +            front :
                  //          / \
                  //         *-- *          *-- *
                  //        /     \
                  //       *       *      *       *
                  //      / \     / \      \     /
                  //     + --*---*-- +      *   *
                  */
                  assert((edgesVLSdata[0]->nbcuts() == 2) && (edgesVLSdata[1]->nbcuts() == 2) &&
                         (edgesVLSdata[2]->nbcuts() == 2));
                  AOMD::mVertex *vT0 = static_cast<AOMD::mVertex *>(face->get(0, 0));
                  AOMD::mVertex *vT1 = static_cast<AOMD::mVertex *>(face->get(0, 1));
                  AOMD::mVertex *vT2 = static_cast<AOMD::mVertex *>(face->get(0, 2));
                  AOMD::mVertex *ve0_0 = (*edgesVLSdata[0])(0);
                  AOMD::mVertex *ve0_1 = (*edgesVLSdata[0])(1);
                  if (vT0 != face->get(1, 0)->get(0, 0)) swap(ve0_0, ve0_1);
                  AOMD::mVertex *ve1_0 = (*edgesVLSdata[1])(0);
                  AOMD::mVertex *ve1_1 = (*edgesVLSdata[1])(1);
                  if (vT1 != face->get(1, 1)->get(0, 0)) swap(ve1_0, ve1_1);
                  AOMD::mVertex *ve2_0 = (*edgesVLSdata[2])(0);
                  AOMD::mVertex *ve2_1 = (*edgesVLSdata[2])(1);
                  if (vT2 != face->get(1, 2)->get(0, 0)) swap(ve2_0, ve2_1);
                  boundarymmesh.createEdge(ve0_1, ve1_0, GEn1);
                  boundarymmesh.createEdge(ve1_1, ve2_0, GEn1);
                  boundarymmesh.createEdge(ve2_1, ve0_0, GEn1);
                  triangle Tout0 = {vT0, ve0_0, ve2_1};
                  triangle Tout1 = {vT1, ve1_0, ve0_1};
                  triangle Tout2 = {vT2, ve2_0, ve1_1};
                  triangle Tin0 = {ve0_0, ve2_0, ve2_1};
                  triangle Tin1 = {ve0_0, ve1_1, ve2_0};
                  triangle Tin2 = {ve0_0, ve1_0, ve1_1};
                  triangle Tin3 = {ve0_0, ve0_1, ve1_0};
                  tcutdata.TrianglesOut = {Tout0, Tout1, Tout2};
                  tcutdata.TrianglesIn = {Tin0, Tin1, Tin2, Tin3};
                  break;
               }
               default:
                  throw;
            }  // end switch case 3
            break;
         }
         default:
            throw;
      }  // end gen swith
}

#ifdef WITH_XLEGACYSIMPLECUT
void xPhysSurfVLS::tagAndPromoteEntity(AOMD::mEntity *e, char tag)
{
   filter_data.setData(*e) = tag;
   for_each(promotors.begin(), promotors.end(), std::bind(&xPhysSurfByTagging::promoteStatus, std::placeholders::_1, e));
}

void xPhysSurfVLS::unTagAndUnPromoteEntity(AOMD::mEntity *e)
{
   filter_data.deleteData(*e);
   for_each(promotors.begin(), promotors.end(), std::bind(&xPhysSurfByTagging::unPromoteStatus, std::placeholders::_1, e));
}

void xPhysSurfVLS::treatCutData(AOMD::mEntity *face, xVLSTriangleCutData &fcutdata)
{
   if (fcutdata.fromMeshIn)
   {
      xfem::xPartition part;
      xfem::xMesh::getPartition(face, part, xfem::xAcceptAll());
      for (AOMD::mEntity *pe : part) filter_data.setData(*pe) = 'i';
      filter_data.setData(*face) = 'i';
   }
   if (fcutdata.fromMeshOut)
   {
      xfem::xPartition part;
      xfem::xMesh::getPartition(face, part, xfem::xAcceptAll());
      for (AOMD::mEntity *pe : part) filter_data.setData(*pe) = 'o';
      filter_data.setData(*face) = 'o';
   }

   const auto sizeIn = fcutdata.TrianglesIn.size();
   const auto sizeOut = fcutdata.TrianglesOut.size();

   if ((sizeIn || sizeOut))
   {
      filter_data.setData(*face) = 'c';
      xVLSTriangleCutAttachableData &fdata = faces_data.setData(*face);
      for_each(promotors.begin(), promotors.end(),
               std::bind(&xPhysSurfByTagging::setPromotorStatus, std::placeholders::_1, face));
      if (sizeIn)
      {
         std::set<AOMD::mEntity *> &entities = fdata.TrianglesIn;
         std::transform(fcutdata.TrianglesIn.begin(), fcutdata.TrianglesIn.end(), inserter(entities, entities.begin()),
                        std::bind(new_mEntityFromTriangle, std::placeholders::_1, face->getClassification()));
         for_each(entities.begin(), entities.end(),
                  std::bind(&xPhysSurfVLS::tagAndPromoteEntity, this, std::placeholders::_1, 'i'));
      }
      if (sizeOut)
      {
         std::set<AOMD::mEntity *> &entities = fdata.TrianglesOut;
         ;
         std::transform(fcutdata.TrianglesOut.begin(), fcutdata.TrianglesOut.end(), inserter(entities, entities.begin()),
                        std::bind(new_mEntityFromTriangle, std::placeholders::_1, face->getClassification()));

         for_each(entities.begin(), entities.end(),
                  std::bind(&xPhysSurfVLS::tagAndPromoteEntity, this, std::placeholders::_1, 'o'));
      }
   }
}

void xPhysSurfVLS::clearCutData(AOMD::mEntity *face)
{
   filter_data.deleteData(*face);
   xfem::xPartition part;
   xfem::xMesh::getPartition(face, part, xfem::xAcceptAll());
   for (auto pe : part) filter_data.deleteData(*pe);
   xVLSTriangleCutAttachableData *cutdata = faces_data.getData(*face);
   if (cutdata)
   {
      for_each(promotors.begin(), promotors.end(),
               std::bind(&xPhysSurfByTagging::setPromotorStatus, std::placeholders::_1, face));
      for_each(cutdata->TrianglesIn.begin(), cutdata->TrianglesIn.end(),
               std::bind(&xPhysSurfVLS::unTagAndUnPromoteEntity, this, std::placeholders::_1));
      for_each(cutdata->TrianglesOut.begin(), cutdata->TrianglesOut.end(),
               std::bind(&xPhysSurfVLS::unTagAndUnPromoteEntity, this, std::placeholders::_1));
      for_each(cutdata->TrianglesIn.begin(), cutdata->TrianglesIn.end(), delete_ptr());
      for_each(cutdata->TrianglesOut.begin(), cutdata->TrianglesOut.end(), delete_ptr());
      faces_data.deleteData(*face);
   }
}

bool xPhysSurfVLS::isIn(const AOMD::mEntity &e) const
{
   const char *ftag = filter_data.getData(e);
   return ftag ? (*ftag == 'i') : false;
}
bool xPhysSurfVLS::isCut(const AOMD::mEntity &e) const
{
   const char *ftag = filter_data.getData(e);
   return ftag ? (*ftag == 'c') : false;
}
bool xPhysSurfVLS::isOut(const AOMD::mEntity &e) const
{
   const char *ftag = filter_data.getData(e);
   return ftag ? (*ftag == 'o') : false;
}

bool xPhysSurfVLS::isSupportCoversIn(AOMD::mEntity *e) const
{
   for (auto esup : getsupport(e))
      if (isIn(*esup) || isCut(*esup)) return true;
   return false;
}

bool xPhysSurfVLS::isSupportCoversOut(AOMD::mEntity *e) const
{
   for (auto esup : getsupport(e))
      if (isOut(*esup) || isCut(*esup)) return true;
   return false;
}

xfem::xEntityFilter xPhysSurfVLS::isSupportCoversOutFilter() const
{
   return bind1st(mem_fun(&xPhysSurfVLS::isSupportCoversOut), this);
}

xfem::xEntityFilter xPhysSurfVLS::isSupportCoversInFilter() const
{
   return bind1st(mem_fun(&xPhysSurfVLS::isSupportCoversIn), this);
}

xfem::xEntityFilter xPhysSurfVLS::isInFilter() const
{
   return [this](AOMD::mEntity *e) { return this->isIn(*e); };
}

xfem::xEntityFilter xPhysSurfVLS::isCutFilter() const
{
   return [this](AOMD::mEntity *e) { return this->isCut(*e); };
}
xfem::xEntityFilter xPhysSurfVLS::isOutFilter() const
{
   return [this](AOMD::mEntity *e) { return this->isOut(*e); };
}

xPhysSurfVLS::~xPhysSurfVLS()
{
   std::for_each(facestoclean.begin(), facestoclean.end(), std::bind(&xPhysSurfVLS::clearCutData, this, std::placeholders::_1));
   clearAllSupportComponent();
}

void xPhysSurfVLS::getPartitionOut(AOMD::mEntity *e, xfem::xPartition &partition, xfem::xEntityFilter filter)
{
   xVLSTriangleCutAttachableData *cutdata = faces_data.getData(*e);
   if (cutdata)
   {
      xVLSTriangleCutAttachableData::entityContainer TrianglesOut = cutdata->TrianglesOut;
      xVLSTriangleCutAttachableData::entityContainer::iterator itTout = TrianglesOut.begin();
      xVLSTriangleCutAttachableData::entityContainer::iterator itToutend = TrianglesOut.end();
      for (; itTout != itToutend; ++itTout)
      {
         if (filter(*itTout)) partition.insert(*itTout);
      }
   }
   else if (isOut(*e))
      xfem::xMesh::getPartition(e, partition, filter);
   //   else xMesh::getPartition( e,  partition, filter); //this is an ugly hack for some boundary conditions ..
}

void xPhysSurfVLS::getPartition(AOMD::mEntity *e, xfem::xPartition &partition, xfem::xEntityFilter filter)
{
   xVLSTriangleCutAttachableData *cutdata = faces_data.getData(*e);
   if (cutdata)
   {
      xVLSTriangleCutAttachableData::entityContainer &TrianglesOut = cutdata->TrianglesOut;
      xVLSTriangleCutAttachableData::entityContainer::iterator itTout = TrianglesOut.begin();
      xVLSTriangleCutAttachableData::entityContainer::iterator itToutend = TrianglesOut.end();
      for (; itTout != itToutend; ++itTout)
      {
         if (filter(*itTout)) partition.insert(*itTout);
      }
      xVLSTriangleCutAttachableData::entityContainer &TrianglesIn = cutdata->TrianglesIn;
      xVLSTriangleCutAttachableData::entityContainer::iterator itTin = TrianglesIn.begin();
      xVLSTriangleCutAttachableData::entityContainer::iterator itTinend = TrianglesIn.end();
      for (; itTin != itTinend; ++itTin)
      {
         if (filter(*itTin)) partition.insert(*itTin);
      }
   }
   else
      xfem::xMesh::getPartition(e, partition, filter);
}

void xPhysSurfVLS::clearSupportComponent(AOMD::mEntity *e)
{
   supportcomponent_data.deleteData(*e);
   knowncomponent.erase(e);
}

void xPhysSurfVLS::clearAllSupportComponent()
{
   for (AOMD::mEntity *e : knowncomponent) clearSupportComponent(e);
}

const xfem::xSupportComponent &xPhysSurfVLS::getSupportComponent(AOMD::mEntity *e) const
{
   xfem::xSupportComponent *psupportcomponents = supportcomponent_data.getData(*e);
   if (!psupportcomponents)
   {
      psupportcomponents = &supportcomponent_data.setData(*e);
      knowncomponent.insert(e);
      std::set<AOMD::mEntity *> support = getsupport(e);
      size_t size_support = support.size();
      if (dim != 2) throw;
      std::vector<const AOMD::mEntity *> TrianglesOut;
      xfem::xSupportComponent::component &deadzone = psupportcomponents->getSplitZone();
      TrianglesOut.reserve(4 * size_support);
      double vol_support = 0.;
      for (const auto ce : support)
      {
         xfem::xElement elce(const_cast<AOMD::mEntity *>(ce));
         vol_support += elce.getVolume();
         if (const xVLSTriangleCutAttachableData *cutdata = faces_data.getData(*ce))
         {
            const xVLSTriangleCutAttachableData::entityContainer &TrianglesOutLoc = cutdata->TrianglesOut;
            TrianglesOut.insert(TrianglesOut.end(), TrianglesOutLoc.begin(), TrianglesOutLoc.end());
            const xVLSTriangleCutAttachableData::entityContainer &TrianglesInLoc = cutdata->TrianglesIn;
            deadzone.insert(deadzone.end(), TrianglesInLoc.begin(), TrianglesInLoc.end());
         }
         else
         {
            if (isOut(*ce)) TrianglesOut.push_back(ce);
            if (isIn(*ce)) deadzone.push_back(ce);
         }
      }
      std::vector<AOMD::mVertex *> grapid_to_vertex;
      datamanager<size_t> vertex_to_id;
      size_t id = 0;
      typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
      Graph G;
      for (const AOMD::mEntity *Triangle : TrianglesOut)
      {
         size_t vids[3];
         for (int i = 0; i < 3; ++i)
         {
            AOMD::mVertex *v = static_cast<AOMD::mVertex *>(Triangle->get(0, i));
            size_t *pvid = vertex_to_id.getData(*v);
            if (!pvid)
            {
               vertex_to_id.setData(*v) = id;
               grapid_to_vertex.push_back(v);
               vids[i] = id;
               ++id;
            }
            else
               vids[i] = *pvid;
         }
         add_edge(vids[0], vids[1], G);
         add_edge(vids[1], vids[2], G);
      }
      std::vector<int> component(num_vertices(G));
      size_t num = size_t(boost::connected_components(G, &component[0]));
      std::vector<xfem::xSupportComponent::component> componenttab(num);
      std::vector<double> component_vol(num, 0.);
      for (const AOMD::mEntity *Triangle : TrianglesOut)
      {
         size_t vid = vertex_to_id.at(*Triangle->get(0, 0));
         size_t componentid = size_t(component[vid]);
         componenttab[componentid].push_back(Triangle);
         xfem::xElement elT(const_cast<AOMD::mEntity *>(Triangle));
         component_vol[componentid] += elT.getVolume();
      }
      double eps_vol_component = eps4 * vol_support;

      for (size_t i = 0; i < num; ++i)
      {
         if (component_vol[i] > eps_vol_component)
            psupportcomponents->addComponent(componenttab[i]);
         else
         {
            deadzone.insert(deadzone.end(), componenttab[i].begin(), componenttab[i].end());
         }
      }
   }
   return *psupportcomponents;
}
#endif

double areaTri(AOMD::mEntity *e)
{
   const xtensor::xPoint A = static_cast<AOMD::mVertex *>(e->get(0, 0))->point();
   const xtensor::xPoint B = static_cast<AOMD::mVertex *>(e->get(0, 1))->point();
   const xtensor::xPoint C = static_cast<AOMD::mVertex *>(e->get(0, 2))->point();
   return 0.5 * (xtensor::xVector<>(A, B) % xtensor::xVector<>(A, C)).mag();
}

}  // end namespace xcut
