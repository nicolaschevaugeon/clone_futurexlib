/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _SCAN_ELEMENT_IMP_H
#define _SCAN_ELEMENT_IMP_H
#ifdef DEBUG
#include <iostream>
#endif

namespace xgeom
{
template <typename UNARYOP>
xScanElement2D<UNARYOP>::xScanElement2D(typename UNARYOP::PointArg_t p_min_, typename UNARYOP::PointArg_t p_max_,
                                        const int level_, const double tolerance_)
    : p_min(entity_to_scan_entity.convertPoint(p_min_)),
      p_max(entity_to_scan_entity.convertPoint(p_max_)),
      level(level_),
      tolerance(tolerance_),
      nb(int(pow(2., level) + 1)),
      pas_x((p_max[0] - p_min[0]) / (nb - 1)),
#ifdef DEBUG
      pas_y((p_max[1] - p_min[1]) / (nb - 1)),
      debug(false)
#else
      pas_y((p_max[1] - p_min[1]) / (nb - 1))
#endif
{
#ifdef DEBUG
   assert(level > 0);
   assert(pas_x > tolerance);
   assert(pas_y > tolerance);
#endif
}

template <typename UNARYOP>
inline typename xScanElement2D<UNARYOP>::IJK xScanElement2D<UNARYOP>::getIJK(const int id)
{
   IJK ijk;
   int j = int(id / nb);
   int i = id - j * nb;
   ijk.push_back(i);
   ijk.push_back(j);
   return ijk;
}

template <typename UNARYOP>
inline typename xScanElement2D<UNARYOP>::XYZ xScanElement2D<UNARYOP>::getXYZ(const int id)
{
   XYZ xyz;
   int j = int(id / nb);
   int i = id - j * nb;
   xyz.push_back(p_min[0] + i * pas_x);
   xyz.push_back(p_min[1] + j * pas_y);
   xyz.push_back(0.);
   return xyz;
}

template <typename UNARYOP>
inline double xScanElement2D<UNARYOP>::intersectionLineLine(const double x1, const double y1, const double x2, const double y2,
                                                            const double y_line)
{
   return ((x1 - x2) * y_line - (y2 * x1 - y1 * x2)) / (y1 - y2);
}

template <typename UNARYOP>
void xScanElement2D<UNARYOP>::scanPoint(const double x, const int j, std::vector<int>& active_id)
{
   double divx = (x - p_min[0]) / pas_x;

   int nx_min = (int)floor(divx);
   int nx_max = (int)ceil(divx);

   double x_min = p_min[0] + nx_min * pas_x;
   double x_max = p_min[0] + nx_max * pas_x;

   if (fabs(x - x_min) < tolerance)
   {
      int id = nx_min + j * nb;
      active_id.push_back(id);
   }
   else if (fabs(x - x_max) < tolerance)
   {
      int id = nx_max + j * nb;
      active_id.push_back(id);
   }
}

template <typename UNARYOP>
void xScanElement2D<UNARYOP>::scan(typename UNARYOP::ElementArg_t triangle_, std::vector<int>& active_id)
{
   Element_t triangle = entity_to_scan_entity.convertElement(triangle_);
   std::stable_sort(triangle.begin(), triangle.end(), xGreaterPoint2D());

   // Points classes
   Point_t p1 = triangle[0];
   Point_t p2 = triangle[1];
   Point_t p3 = triangle[2];

   double y_max = p1[1];
   double y_min = p3[1];

   double div1 = (y_max - p_min[1]) / pas_y;
   double div2 = (y_min - p_min[1]) / pas_y;

   int n_max = (int)ceil(div1);
   int n_min = (int)floor(div2);  // n_min <= y_min < y_max <= n_max

#ifdef DEBUG
   assert(n_min < n_max);  // must have this
#endif

   for (int j = n_max; j >= n_min; --j)  // j>n_min has been changed => it seems to correct a bug on y_line = 0
   {
      double y_line = p_min[1] + j * pas_y;

      if ((y_line < (y_max + tolerance)) && (y_line > (y_min - tolerance)))
      {
         double x1 = p1[0];
         double y1 = p1[1];
         double x2 = p2[0];
         double y2 = p2[1];
         double x3 = p3[0];
         double y3 = p3[1];

         double xR, xL;
         bool compute = true;

         if ((fabs(y_line - y1) < tolerance) && (fabs(y_line - y2) < tolerance))  // (1) highest and middle nodes on line
         {
#ifdef DEBUG
            if (debug) std::cout << "                                                 < 1, 2 on line > " << std::endl;
#endif
            xL = x1;
            xR = x2;
         }
         else if ((fabs(y_line - y2) < tolerance) && (fabs(y_line - y3) < tolerance))  // (2) middle and lowest nodes on line
         {
#ifdef DEBUG
            if (debug) std::cout << "                                                 < 2, 3 on line > " << std::endl;
#endif
            xL = x2;
            xR = x3;
         }
         else if (fabs(y_line - y1) < tolerance)  // (3) highest node on line only
         {
#ifdef DEBUG
            if (debug) std::cout << "                                                 < 1 on line NO COMPUTATION > " << std::endl;
#endif
            compute = false;
            scanPoint(x1, j, active_id);
         }
         else if (fabs(y_line - y3) < tolerance)  // (4) lowest node on line only
         {
#ifdef DEBUG
            if (debug) std::cout << "                                                 < 3 on line NO COMPUTATION > " << std::endl;
#endif
            compute = false;
            scanPoint(x3, j, active_id);
         }
         else if (fabs(y_line - y2) < tolerance)  // (5) middle node on line
         {
#ifdef DEBUG
            if (debug) std::cout << "                                                 < 2 on line > " << std::endl;
#endif
            xL = x2;
            xR = intersectionLineLine(x1, y1, x3, y3, y_line);
         }
         else if (y_line > y2)  // (6) no node on edge but line over middle node
         {
#ifdef DEBUG
            if (debug) std::cout << "                                                 < 2 under line > " << std::endl;
#endif
            xL = intersectionLineLine(x1, y1, x2, y2, y_line);
            xR = intersectionLineLine(x1, y1, x3, y3, y_line);
         }
         else if (y_line < y2)  // (7) no node on edge but line under middle node
         {
#ifdef DEBUG
            if (debug) std::cout << "                                                 < 2 over line > " << std::endl;
#endif
            xL = intersectionLineLine(x1, y1, x3, y3, y_line);
            xR = intersectionLineLine(x2, y2, x3, y3, y_line);
         }
         else
         {
            std::cerr << "Bug in xScanElement_imp.h" << std::endl;
            throw;  // Other situation ? Bug !!!
         }

         if (compute)
         {
            if (xL > xR)
            {
               double buf = xL;
               xL = xR;
               xR = buf;
            }
#ifdef DEBUG
            if (debug) std::cout << "x_min " << xL << " x_max " << xR << std::endl;
            assert(xL < xR);
#endif
            double div3 = (xR - p_min[0]) / pas_x;
            double div4 = (xL - p_min[0]) / pas_x;

            int nx_max = (int)ceil(div3);
            int nx_min = (int)floor(div4);

            if ((xL - nx_min * pas_x) > tolerance)  // on reduit a gauche
            {
               nx_min++;
            }
            if ((nx_max * pas_x - xR) > tolerance)  // on reduit a droite
            {
               nx_max--;
            }

            // if (nx_min <= nx_max) { // securise

            for (int i = nx_min; i <= nx_max; ++i)
            {
               int id = i + j * nb;
               active_id.push_back(id);
            }
            //}
         }
      }
   }
}

template <typename UNARYOP>
xScanElement3D<UNARYOP>::xScanElement3D(typename UNARYOP::PointArg_t p_min_, typename UNARYOP::PointArg_t p_max_,
                                        const int level_, const double tolerance_)
    : p_min(entity_to_scan_entity.convertPoint(p_min_)),
      p_max(entity_to_scan_entity.convertPoint(p_max_)),
      level(level_),
      tolerance(tolerance_),
      nb(int(pow(2., level) + 1)),
      nb2(nb * nb),
      pas_x((p_max[0] - p_min[0]) / (nb - 1)),
      pas_y((p_max[1] - p_min[1]) / (nb - 1)),
      pas_z((p_max[2] - p_min[2]) / (nb - 1)),
#ifdef DEBUG
      debug(false),
#endif
      scan2D(p_min, p_max, level, tolerance)
{
#ifdef DEBUG
   assert(level > 0);
   assert(pas_x > tolerance);
   assert(pas_y > tolerance);
   assert(pas_z > tolerance);
#endif
}

template <typename UNARYOP>
inline typename xScanElement3D<UNARYOP>::IJK xScanElement3D<UNARYOP>::getIJK(const int id)
{
   IJK ijk;
   int k = int(id / nb2);
   int rest = id - k * nb2;
   int j = int(rest / nb);
   int i = id - j * nb - k * nb2;
   ijk.push_back(i);
   ijk.push_back(j);
   ijk.push_back(k);
   return ijk;
}

template <typename UNARYOP>
inline typename xScanElement3D<UNARYOP>::XYZ xScanElement3D<UNARYOP>::getXYZ(const int id)
{
   XYZ xyz;
   int k = int(id / nb2);
   int rest = id - k * nb2;
   int j = int(rest / nb);
   int i = id - j * nb - k * nb2;
   xyz.push_back(p_min[0] + i * pas_x);
   xyz.push_back(p_min[1] + j * pas_y);
   xyz.push_back(p_min[2] + k * pas_z);
   return xyz;
}

template <typename UNARYOP>
inline typename UNARYOP::PointRes_t xScanElement3D<UNARYOP>::intersectionLinePlane(const double x1, const double y1,
                                                                                   const double z1, const double x2,
                                                                                   const double y2, const double z2,
                                                                                   const double z_plane)
{
   double ratio = (z1 - z_plane) / (z1 - z2);
   double x = x1 + (x2 - x1) * ratio;
   double y = y1 + (y2 - y1) * ratio;
   typename UNARYOP::PointRes_t point;
   point.push_back(x);
   point.push_back(y);
   return point;
}

template <typename UNARYOP>
void xScanElement3D<UNARYOP>::scanPoint(const double x, const double y, const int k, std::vector<int>& active_id)
{
   double divx = (x - p_min[0]) / pas_x;
   double divy = (y - p_min[1]) / pas_y;

   int nx_min = (int)floor(divx);
   int nx_max = (int)ceil(divx);
   int ny_min = (int)floor(divy);
   int ny_max = (int)ceil(divy);

   double x_min = p_min[0] + nx_min * pas_x;
   double x_max = p_min[0] + nx_max * pas_x;
   double y_min = p_min[1] + ny_min * pas_y;
   double y_max = p_min[1] + ny_max * pas_y;

   if (fabs(x - x_min) < tolerance)
   {
      if (fabs(y - y_min) < tolerance)
      {
         int id = nx_min + ny_min * nb + k * nb2;
         active_id.push_back(id);
      }
      else if (fabs(y - y_max) < tolerance)
      {
         int id = nx_min + ny_max * nb + k * nb2;
         active_id.push_back(id);
      }
   }
   else if (fabs(x - x_max) < tolerance)
   {
      if (fabs(y - y_min) < tolerance)
      {
         int id = nx_max + ny_min * nb + k * nb2;
         active_id.push_back(id);
      }
      else if (fabs(y - y_max) < tolerance)
      {
         int id = nx_max + ny_max * nb + k * nb2;
         active_id.push_back(id);
      }
   }
}

template <typename UNARYOP>
void xScanElement3D<UNARYOP>::scanSegment(const double x1, const double y1, const double x2, const double y2, const int k,
                                          std::vector<int>& active_id)
{
   double y_max = std::max(y1, y2);
   double y_min = std::min(y1, y2);

   double div1 = (y_max - p_min[1]) / pas_y;
   double div2 = (y_min - p_min[1]) / pas_y;

   int ny_max = (int)ceil(div1);
   int ny_min = (int)floor(div2);

   if (fabs(y_max - y_min) < tolerance)  // (1) Segment parallel to y
   {
      double x_max = std::max(x1, x2);
      double x_min = std::min(x1, x2);

      double div1x = (x_max - p_min[0]) / pas_x;
      double div2x = (x_min - p_min[0]) / pas_x;

      int nx_max = (int)ceil(div1x);
      int nx_min = (int)floor(div2x);

#ifdef DEBUG
      assert(nx_min < nx_max);
#endif

      if ((x_min - (p_min[0] + nx_min * pas_x)) > tolerance)  // on reduit a gauche
      {
         nx_min++;
      }
      if (((p_min[0] + nx_max * pas_x) - x_max) > tolerance)  // on reduit a droite
      {
         nx_max--;
      }

      for (int i = nx_min; i <= nx_max; ++i)
      {
         int id = i + ny_max * nb + k * nb2;
         active_id.push_back(id);
      }
   }
   else  // (2) Otherwise
   {
#ifdef DEBUG
      assert(ny_min < ny_max);  // must be verify
#endif

      for (int j = ny_max; j > ny_min; --j)
      {
         double y_line = p_min[1] + j * pas_y;

         if ((y_line < (y_max + tolerance)) && (y_line > (y_min - tolerance)))
         {
            double x = scan2D.intersectionLineLine(x1, y1, x2, y2, y_line);
            std::vector<int> active_id_2D;
            scan2D.scanPoint(x, j, active_id_2D);

            for (std::vector<int>::iterator it = active_id_2D.begin(); it != active_id_2D.end(); ++it)
            {
               int id = *it;
               id += k * nb2;
               active_id.push_back(id);
            }
         }
      }
   }
}

template <typename UNARYOP>
void xScanElement3D<UNARYOP>::scan(typename UNARYOP::ElementArg_t tetraedre_, std::vector<int>& active_id)
{
   Element_t tetraedre = entity_to_scan_entity.convertElement(tetraedre_);
   std::stable_sort(tetraedre.begin(), tetraedre.end(), xGreaterPoint3D());

   // Points classes
   Point_t p1 = tetraedre[0];
   Point_t p2 = tetraedre[1];
   Point_t p3 = tetraedre[2];
   Point_t p4 = tetraedre[3];
#ifdef DEBUG
   if (debug) std::cout << "--> in scanConvert" << std::endl;
#endif
   double z_max = p1[2];
   double z_min = p4[2];
#ifdef DEBUG
   if (debug) std::cout << "z_min " << z_min << " z_max " << z_max << std::endl;
#endif
   double div1 = (z_max - p_min[2]) / pas_z;
   double div2 = (z_min - p_min[2]) / pas_z;

   int n_max = (int)ceil(div1);
   int n_min = (int)floor(div2);
#ifdef DEBUG
   if (debug) std::cout << "n_min " << n_min << " n_max " << n_max << std::endl;
      // assert(n_min < n_max);
#endif
   for (int k = n_max; k > n_min; --k)
   {
      double z_plane = p_min[2] + k * pas_z;

      if ((z_plane < (z_max + tolerance)) && (z_plane > (z_min - tolerance)))
      {
         double x1 = p1[0];
         double y1 = p1[1];
         double z1 = p1[2];
         double x2 = p2[0];
         double y2 = p2[1];
         double z2 = p2[2];
         double x3 = p3[0];
         double y3 = p3[1];
         double z3 = p3[2];
         double x4 = p4[0];
         double y4 = p4[1];
         double z4 = p4[2];

         bool compute = true;
         Element_t triangle;
         Element_t triangle2;

         if ((fabs(z_plane - z1) < tolerance) && (fabs(z_plane - z2) < tolerance) &&
             (fabs(z_plane - z3) < tolerance))  // (1) 1, 2 and 3 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 1, 2 and 3 on plane > " << std::endl;
#endif
            Point_t pt1;
            Point_t pt2;
            Point_t pt3;
            pt1.push_back(x1);
            pt1.push_back(y1);
            pt2.push_back(x2);
            pt2.push_back(y2);
            pt3.push_back(x3);
            pt3.push_back(y3);
            triangle.push_back(pt1);
            triangle.push_back(pt2);
            triangle.push_back(pt3);
         }
         else if ((fabs(z_plane - z2) < tolerance) && (fabs(z_plane - z3) < tolerance) &&
                  (fabs(z_plane - z4) < tolerance))  // (2) 2, 3 and 4 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 2, 3 and 4 on plane > " << std::endl;
#endif
            Point_t pt2;
            Point_t pt3;
            Point_t pt4;
            pt2.push_back(x2);
            pt2.push_back(y2);
            pt3.push_back(x3);
            pt3.push_back(y3);
            pt4.push_back(x4);
            pt4.push_back(y4);
            triangle.push_back(pt2);
            triangle.push_back(pt3);
            triangle.push_back(pt4);
         }
         else if ((fabs(z_plane - z1) < tolerance) && (fabs(z_plane - z2) < tolerance))  // (3) 1 and 2 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 1, 2 on plane NO COMPUTING > " << std::endl;
#endif
            compute = false;
            scanSegment(x1, y1, x2, y2, k, active_id);
         }
         else if ((fabs(z_plane - z3) < tolerance) && (fabs(z_plane - z4) < tolerance))  // (4) 3 and 4 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 3, 4 on plane NO COMPUTING > " << std::endl;
#endif
            compute = false;
            scanSegment(x3, y3, x4, y4, k, active_id);
         }
         else if ((fabs(z_plane - z1) < tolerance))  // (5) 1 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 1 on plane NO COMPUTING > " << std::endl;
#endif
            compute = false;
            scanPoint(x1, y1, k, active_id);
         }
         else if ((fabs(z_plane - z4) < tolerance))  // (6) 4 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 4 on plane NO COMPUTING > " << std::endl;
#endif
            compute = false;
            scanPoint(x4, y4, k, active_id);
         }
         else if ((fabs(z_plane - z2) < tolerance) && (fabs(z_plane - z3) < tolerance))  // (7) 2 and 3 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 2, 3 on plane > " << std::endl;
#endif
            Point_t pt2;
            Point_t pt3;
            pt2.push_back(x2);
            pt2.push_back(y2);
            pt3.push_back(x3);
            pt3.push_back(y3);
            Point_t pt_4 = intersectionLinePlane(x1, y1, z1, x4, y4, z4, z_plane);  // 1-x-4
            triangle.push_back(pt2);
            triangle.push_back(pt3);
            triangle.push_back(pt_4);
         }
         else if ((fabs(z_plane - z2) < tolerance))  // (8) 2 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 2 on plane > " << std::endl;
#endif
            Point_t pt2;
            pt2.push_back(x2);
            pt2.push_back(y2);
            Point_t pt_1 = intersectionLinePlane(x1, y1, z1, x3, y3, z3, z_plane);  // 1-x-3
            Point_t pt_2 = intersectionLinePlane(x1, y1, z1, x4, y4, z4, z_plane);  // 1-x-4
            triangle.push_back(pt2);
            triangle.push_back(pt_1);
            triangle.push_back(pt_2);
         }
         else if ((fabs(z_plane - z3) < tolerance))  // (9) 3 on plane
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < 3 on plane > " << std::endl;
#endif
            Point_t pt3;
            pt3.push_back(x3);
            pt3.push_back(y3);
            Point_t pt_1 = intersectionLinePlane(x1, y1, z1, x4, y4, z4, z_plane);  // 1-x-4
            Point_t pt_2 = intersectionLinePlane(x2, y2, z2, x4, y4, z4, z_plane);  // 2-x-4
            triangle.push_back(pt3);
            triangle.push_back(pt_1);
            triangle.push_back(pt_2);
         }
         else if (z_plane > z2)  // (10) no node on plane but plane over 2
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < plane over 2 > " << std::endl;
#endif
            Point_t pt_1 = intersectionLinePlane(x1, y1, z1, x2, y2, z2, z_plane);  // 1-x-2
            Point_t pt_2 = intersectionLinePlane(x1, y1, z1, x3, y3, z3, z_plane);  // 1-x-3
            Point_t pt_3 = intersectionLinePlane(x1, y1, z1, x4, y4, z4, z_plane);  // 1-x-4
            triangle.push_back(pt_1);
            triangle.push_back(pt_2);
            triangle.push_back(pt_3);
         }
         else if (z_plane < z3)  // (11) no node on plane but plane under 3
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < plane under 3 > " << std::endl;
#endif
            Point_t pt_1 = intersectionLinePlane(x1, y1, z1, x4, y4, z4, z_plane);  // 1-x-4
            Point_t pt_2 = intersectionLinePlane(x2, y2, z2, x4, y4, z4, z_plane);  // 2-x-4
            Point_t pt_3 = intersectionLinePlane(x3, y3, z3, x4, y4, z4, z_plane);  // 3-x-4
            triangle.push_back(pt_1);
            triangle.push_back(pt_2);
            triangle.push_back(pt_3);
         }
         else if ((z_plane < z2) &&
                  (z_plane > z3))  // (12) no node on plane but plane between 2 and 3 : careful two triangles !!!
         {
#ifdef DEBUG
            if (debug) std::cout << "                                       < plane between 2, 3 > " << std::endl;
#endif
            Point_t pt_1 = intersectionLinePlane(x1, y1, z1, x3, y3, z3, z_plane);  // 1-x-3
            Point_t pt_2 = intersectionLinePlane(x1, y1, z1, x4, y4, z4, z_plane);  // 1-x-4
            Point_t pt_3 = intersectionLinePlane(x2, y2, z2, x3, y3, z3, z_plane);  // 2-x-3
            Point_t pt_4 = intersectionLinePlane(x2, y2, z2, x4, y4, z4, z_plane);  // 2-x-4
            triangle.push_back(pt_2);
            triangle.push_back(pt_4);
            triangle.push_back(pt_3);
            triangle2.push_back(pt_1);
            triangle2.push_back(pt_2);
            triangle2.push_back(pt_4);
         }
         else
         {
            std::cerr << "Bug in xScanElement_imp.h" << std::endl;
            throw;  // Other situation ? Bug !!!
         }

         if (compute)
         {
            std::vector<int> active_id_2D;
            scan2D.scan(triangle, active_id_2D);

            if (!triangle2.empty())
            {
               scan2D.scan(triangle2, active_id_2D);
            }

            for (std::vector<int>::iterator it = active_id_2D.begin(); it != active_id_2D.end(); ++it)
            {
               int id = *it;
               id += k * nb2;
               active_id.push_back(id);
            }
         }
      }
   }
}
}  // namespace xgeom

#endif
