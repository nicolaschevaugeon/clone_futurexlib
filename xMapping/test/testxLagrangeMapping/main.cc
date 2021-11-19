#include <iomanip>
#include <iostream>

// xmapping
#include "xLagrangeMapping.h"

const double eps = 1.e-6;
using point = xtensor::xPoint;
using vect = xmapping::vector3d;
using tens = xmapping::tensor;

std::ostream &operator<<(std::ostream &out, const vect &v)
{
   out << "{" << v(0) << ", " << v(1) << ", " << v(2) << "}";
   return out;
}

std::ostream &operator<<(std::ostream &out, const tens &T)
{
   out << "{{" << T(0, 0) << ", " << T(0, 1) << ", " << T(0, 2) << "},";
   out << "{" << T(1, 0) << ", " << T(1, 1) << ", " << T(1, 2) << "},";
   out << "{" << T(2, 0) << ", " << T(2, 1) << ", " << T(2, 2) << "}}";
   return out;
}

double norm(const tens &in)
{
   return sqrt(in(0, 0) * in(0, 0) + in(0, 1) * in(0, 1) + in(0, 2) * in(0, 2) + in(1, 0) * in(1, 0) + in(1, 1) * in(1, 1) +
               in(1, 2) * in(2, 2) + in(2, 0) * in(2, 0) + in(2, 1) * in(2, 1) + in(2, 2) * in(2, 2));
}

double norm(const vect &in) { return sqrt(in(0) * in(0) + in(1) * in(1) + in(2) * in(2)); }

double norm(const point &in) { return sqrt(in(0) * in(0) + in(1) * in(1) + in(2) * in(2)); }

int main()
{
   std::cout << "Testing xLagrangeMapping" << std::endl;
   bool failed = false;
   // testing edge case.
   {
      std::cout << "Testing for EDGE" << std::endl;
      std::vector<point> p = {{0., 0., 0.}, {1., 2., 3.}};
      xmapping::xLagrangeMapping mapping(xmapping::xReferenceElementType::EDGE, p);
      double u, v, w, x, y, z;
      if (mapping.inReferenceElement(-2., 0., 0.))
      {
         std::cout << "Error evaluating inReferenceElement for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (!mapping.inReferenceElement(0.5, 1., 0.))
      {
         std::cout << "Error evaluating inReferenceElement for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      mapping.COG(u, v, w);
      if (u != 0. || v != 0. || w != 0.)
      {
         std::cout << "Error evaluating COG for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      mapping.eval(u, v, w, x, y, z);
      if (x != 0.5 || y != 1. || z != 1.5)
      {
         std::cout << "Error evaluating eval for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      point min, max;
      point min_ref{0., 0., 0.}, max_ref{1., 2., 3.};
      mapping.boundingBox(min, max);
      if (norm(min - min_ref) > eps || norm(max - max_ref) > eps)
      {
         std::cout << "Error evaluating boudingBox for EDGE" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      const double J_ref = sqrt(14.) / 2.;
      double J_res = mapping.detJac(u, v, w);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating detJac for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens jac;
      J_res = mapping.jacInverse(u, v, w, jac);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens jac_ref{{0.1428571429, -0.894427191, -0.3585685828},
                   {0.2857142857, 0.4472135955, -0.7171371656},
                   {0.4285714286, -0, 0.5976143047}};
      if (norm(jac - jac_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      vect vec_to_push{0.5, 1., 0.};
      vect vec_pushed_res = vec_to_push;
      mapping.pushBack(0., 0., 0., 1, &vec_pushed_res);
      vect vec_pushed_ref = {-0.8229986196, 0.5900707384, 0.2142857143};
      if (norm(vec_pushed_ref - vec_pushed_res) > eps)
      {
         std::cout << "Error evaluating PushBack for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens tens_to_push{{0.5, 1., 0.}, {3., 2., -0.3}, {-1., 15., 7.}};
      tens tens_pushed_res = tens_to_push;
      mapping.pushBack(0., 0., 0., 1, &tens_pushed_res);
      tens tens_pushed_ref{{-2.253284419, -7.024525981, -2.241651922},
                           {2.201635095, -9.576916007, -5.154124238},
                           {-0.3833285904, 9.392785999, 4.183300133}};
      if (norm(tens_pushed_ref - tens_pushed_res) > eps)
      {
         std::cout << "Error evaluating PushBack for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      double uref = 0.2;
      double vref = 0.;
      double wref = 0.;
      u = uref;
      v = vref;
      w = wref;
      mapping.eval(u, v, w, x, y, z);
      mapping.invert(x, y, z, u, v, w);
      if (fabs(uref - u) > eps || fabs(vref - v) > eps || fabs(wref - w) > eps)
      {
         std::cout << "Error evaluating invert for EDGE " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
   }
   // testing triangle case.
   {
      std::cout << "Testing for TRI" << std::endl;
      std::vector<point> p{{0., 0., 0.}, {1., 2., 3.}, {1., 0., 4.}};

      xmapping::xLagrangeMapping mapping(xmapping::xReferenceElementType::TRI, p);
      double u, v, w, x, y, z;
      if (mapping.inReferenceElement(-2., 1., 0.))
      {
         std::cout << "Error evaluating inReferenceElement for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (!mapping.inReferenceElement(0., 0., 0.))
      {
         std::cout << "Error evaluating inReferenceElement for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      if (!mapping.inReferenceElement(0.3, 0.3, 0.))
      {
         std::cout << "Error evaluating inReferenceElement for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      mapping.COG(u, v, w);
      if (fabs(u - 1. / 3.) > eps || fabs(v - 1. / 3.) > eps || fabs(w) > eps)
      {
         std::cout << "Error evaluating COG for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      mapping.eval(u, v, w, x, y, z);
      if (fabs(x - 2. / 3.) > eps || fabs(y - 2. / 3) > eps || fabs(z - 7. / 3.) > eps)
      {
         std::cout << "Error evaluating eval for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      point min, max;
      point min_ref{0., 0., 0.}, max_ref{1., 2., 4.};
      mapping.boundingBox(min, max);
      if (norm(min - min_ref) > eps || norm(max - max_ref) > eps)
      {
         std::cout << "Error evaluating boudingBox for TRI" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      const double J_ref = 8.306623863;
      double J_res = mapping.detJac(u, v, w);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating detJac for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens jac_ref{{0.05797101449, 0.01449275362, 0.9630868247},
                   {0.4927536232, -0.3768115942, -0.1203858531},
                   {-0.01449275362, 0.2463768116, -0.2407717062}};
      tens jac_res;
      J_res = mapping.jacInverse(u, v, w, jac_res);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (norm(jac_res - jac_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      vect vec_to_push{0.5, 1., 0.3};
      vect vec_pushed_ref{0.3324043083, -0.1665505385, 0.1668989229};
      vect vec_pushed_res = vec_to_push;
      mapping.pushBack(0.3, 0.3, 0., 1, &vec_pushed_res);
      if (norm(vec_pushed_res - vec_pushed_ref) > eps)
      {
         std::cout << "Error evaluating PushBack for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens tens_to_push{{0.5, 1., 0.}, {3., 2., -0.3}, {-1., 15., 7.}};
      tens tens_pushed_ref{{-0.8906230566, 14.53325889, 6.737259947},
                           {-0.7636721179, -2.066657362, -0.7296574933},
                           {0.9726557641, -3.133314723, -1.759314987}};
      tens tens_pushed_res = tens_to_push;
      mapping.pushBack(0., 0., 0., 1, &tens_pushed_res);
      if (norm(tens_pushed_res - tens_pushed_ref) > eps)
      {
         std::cout << "Error evaluating PushBack for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      double uref = 0.2;
      double vref = 0.5;
      double wref = 0.;
      u = uref;
      v = vref;
      w = wref;
      mapping.eval(u, v, w, x, y, z);
      mapping.invert(x, y, z, u, v, w);
      if (fabs(uref - u) > eps || fabs(vref - v) > eps || fabs(wref - w) > eps)
      {
         std::cout << "Error evaluating invert for TRI " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }

      {
         vect na;
         // mapping.normalVector(&se10, 0.5,  0., 0. , na);
         const vect ref_n = {-0.0321744726, 0.8365362877, -0.5469660343};
         unsigned short ev[2] = {1, 0};
         vect n = mapping.normalVector(ev, 0.5, 0.);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal TRI " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }

      {
         vect na;
         // mapping.normalVector(&se02, 0.5,  0., 0. , na);
         const vect ref_n = {-0.1167914325, -0.9927271762, 0.02919785812};
         unsigned short ev[2] = {0, 2};
         vect n = mapping.normalVector(ev, 0.5, 0.);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal TRI " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }

      {
         vect na;
         // mapping.normalVector(&se21, 0.5,  0., 0. , na);
         const vect ref_n = {0.269190951, 0.4307055216, 0.8614110433};
         unsigned short ev[2] = {2, 1};
         vect n = mapping.normalVector(ev, 0.5, 0.);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal TRI " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
   }
   // testing quad case.
   {
      std::cout << "Testing for QUAD" << std::endl;
      std::vector<point> p{{0., 0., 0.}, {1., 2., 3.}, {1., 0., 4.}, {0., 0., 4.}};

      /*
     point p0{0., 0., 0.};
     point p1{1., 0., 0.};
     point p2{1., 1., 0.};
     point p3{0., 1., 0.};
     */

      xmapping::xLagrangeMapping mapping(xmapping::xReferenceElementType::QUAD, p);
      double u, v, w, x, y, z;
      if (mapping.inReferenceElement(-2., 1., 0.))
      {
         std::cout << "Error evaluating inReferenceElement for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (!mapping.inReferenceElement(-0.1, -0.2, 0.))
      {
         std::cout << "Error evaluating inReferenceElement for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      if (!mapping.inReferenceElement(0.3, 0.3, 0.))
      {
         std::cout << "Error evaluating inReferenceElement for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      mapping.COG(u, v, w);
      if (u != 0. || v != 0.)
      {
         std::cout << "Error evaluating COG for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      mapping.eval(u, v, w, x, y, z);
      if (fabs(x - 0.5) > eps || fabs(y - 0.5) > eps || fabs(z - 2.75) > eps)
      {
         std::cout << "Error evaluating eval for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      point min, max;
      point min_ref{0., 0., 0.}, max_ref{1., 2., 4.};
      mapping.boundingBox(min, max);
      if (norm(min - min_ref) > eps || norm(max - max_ref) > eps)
      {
         std::cout << "Error evaluating boudingBox for QUAD" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      const double J_ref = 1.205456345;
      double J_res = mapping.detJac(u, v, w);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating detJac for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens jac_ref{{0.623655914, -0.2365591398, 0.8295613558},
                   {0.8602150538, -0.6021505376, -0.5184758474},
                   {0.3440860215, 0.5591397849, -0.2073903389}};
      tens jac_res;
      J_res = mapping.jacInverse(u, v, w, jac_res);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for QUAD" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (norm(jac_res - jac_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      vect vec_to_push{0.5, 1., 0.3};
      vect vec_pushed_res = vec_to_push;
      vect vec_pushed_ref{0.4748094364, -0.3654854803, 0.7438384759};
      mapping.pushBack(0.3, 0.3, 0., 1, &vec_pushed_res);
      if (norm(vec_pushed_res - vec_pushed_ref) > eps)
      {
         std::cout << "Error evaluating PushBack for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens tens_to_push{{0.5, 1., 0.}, {3., 2., -0.3}, {-1., 15., 7.}};
      tens tens_pushed_ref{{-1.227410818, 12.59395797, 5.877897232},
                           {-0.8578682387, -8.121223732, -3.44868577},
                           {2.056852705, -1.648489493, -1.619474308}};
      tens tens_pushed_res = tens_to_push;
      mapping.pushBack(0., 0., 0., 1, &tens_pushed_res);
      if (norm(tens_pushed_res - tens_pushed_ref) > eps)
      {
         std::cout << "Error evaluating PushBack for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      double uref = 0.2;
      double vref = -0.3;
      double wref = 0.;
      u = uref;
      v = vref;
      w = wref;
      mapping.eval(u, v, w, x, y, z);
      mapping.invert(x, y, z, u, v, w);
      if (fabs(uref - u) > eps || fabs(vref - v) > eps)
      {
         std::cout << "Error evaluating invert for QUAD " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }

      {
         vect na;
         // mapping.normalVector(&se10, 0.,  -1. ,0., na);
         const vect ref_n = {0.1741430864, 0.7915594836, -0.5857540178};
         unsigned short ev[2] = {1, 0};
         const vect n = mapping.normalVector(ev, 0., 0.);
         vect e(p[1], p[0]);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na  << " "<< na*e << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << " " << n * e << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal QUAD " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }

      {
         vect na;
         // mapping.normalVector(&se03, -1.,  0., 0. , na);
         const vect ref_n = {-0.7071067812, -0.7071067812, 0};
         unsigned short ev[2] = {0, 3};
         vect n = mapping.normalVector(ev, 0., 0.);
         vect e(p[0], p[3]);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na  << " "<< na*e << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << " " << n * e << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal QUAD " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }

      {
         vect na;
         // mapping.normalVector(&se21, 1.,  0., 0. , na);
         const vect ref_n = {0.4879500365, 0.3903600292, 0.7807200584};
         unsigned short ev[2] = {2, 1};
         const vect n = mapping.normalVector(ev, 0., 0.);
         const vect e(p[2], p[1]);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na  << " "<< na*e << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << " " << n * e << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal QUAD " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
      {
         vect na;
         // mapping.normalVector(&se32, 0.,  1., 0. , na);
         const vect ref_n = {0, -0.3713906764, 0.9284766909};
         //{0.269190951, 0.4307055216, 0.8614110433}
         unsigned short ev[2] = {3, 2};
         vect n = mapping.normalVector(ev, 0., 0.);
         vect e(p[3], p[2]);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na  << " "<< na*e << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << " " << n * e << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal QUAD " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
   }
   // testing tet case.
   {
      std::cout << "Testing for TET" << std::endl;
      std::vector<point> p{{0., 0., 0.}, {1., 0.1, -0.3}, {0.1, 1., -0.1}, {1.1, 0.9, 0.1}};
      xmapping::xLagrangeMapping mapping(xmapping::xReferenceElementType::TET, p);
      double u, v, w, x, y, z;
      if (mapping.inReferenceElement(-2., 1., 0.))
      {
         std::cout << "Error evaluating inReferenceElement for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (mapping.inReferenceElement(-0.1, -0.2, 0.))
      {
         std::cout << "Error evaluating inReferenceElement for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      if (!mapping.inReferenceElement(0.3, 0.3, 0.1))
      {
         std::cout << "Error evaluating inReferenceElement for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      mapping.COG(u, v, w);
      if (u != 0.25 || v != 0.25 || w != 0.25)
      {
         std::cout << "Error evaluating COG for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      mapping.eval(u, v, w, x, y, z);
      if (fabs(x - 0.55) > eps || fabs(y - 0.5) > eps || fabs(z + 0.075) > eps)
      {
         std::cout << "Error evaluating eval for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      point min, max;
      point min_ref{0., 0., -0.3}, max_ref{1.1, 1., 0.1};
      mapping.boundingBox(min, max);
      if (norm(min - min_ref) > eps || norm(max - max_ref) > eps)
      {
         std::cout << "Error evaluating boudingBox for TET" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      const double J_ref = 0.481;
      double J_res = mapping.detJac(u, v, w);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating detJac for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens jac_ref{{0.395010395, -0.5821205821, 0.6029106029},
                   {-0.2494802495, 0.893970894, 0.1455301455},
                   {-2.0997921, -1.642411642, 2.058212058}};
      tens jac_res;
      J_res = mapping.jacInverse(u, v, w, jac_res);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for TET" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (norm(jac_res - jac_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      vect vec_to_push{0.5, 1., 0.3};
      vect vec_pushed_res = vec_to_push;
      vect vec_pushed_ref{-0.2037422037, 0.8128898129, -2.074844075};
      mapping.pushBack(0.3, 0.3, 0., 1, &vec_pushed_res);
      if (norm(vec_pushed_res - vec_pushed_ref) > eps)
      {
         std::cout << "Error evaluating pushBack for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens tens_to_push{{0.5, 1., 0.}, {3., 2., -0.3}, {-1., 15., 7.}};
      tens tens_pushed_ref{{-2.151767152, 8.274428274, 4.395010395},
                           {2.411642412, 3.721413721, 0.7505197505},
                           {-8.035343035, 25.48856549, 14.9002079}};
      tens tens_pushed_res = tens_to_push;
      mapping.pushBack(0., 0., 0., 1, &tens_pushed_res);
      if (norm(tens_pushed_res - tens_pushed_ref) > eps)
      {
         std::cout << "Error evaluating pushBack for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      double uref = 0.2;
      double vref = -0.3;
      double wref = 0.;
      u = uref;
      v = vref;
      w = wref;
      mapping.eval(u, v, w, x, y, z);
      mapping.invert(x, y, z, u, v, w);
      if (fabs(uref - u) > eps || fabs(vref - v) > eps)
      {
         std::cout << "Error evaluating invert for TET " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      {
         std::vector<point> fp = {p[1], p[0], p[2]};
         xmapping::xLagrangeMapping fmapping(xmapping::xReferenceElementType::TRI, fp);
         vect na = fmapping.normalVector(0., 0., 0.);
         const vect ref_n = {-0.280471562, -0.0677000322, -0.957471884};
         unsigned short fv[3] = {1, 0, 2};
         vect nb = mapping.normalVector(fv, 0., 0.);
         std::cout << "Normal n  " << std::setprecision(10) << na << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << nb << std::endl;
         if ((norm(ref_n - na) > eps) || (norm(ref_n - nb) > eps))
         {
            std::cout << "Error evaluating normal TET " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
      {
         vect na;
         std::vector<point> fp = {p[1], p[3], p[2]};
         xmapping::xLagrangeMapping fmapping(xmapping::xReferenceElementType::TRI, fp);
         const vect ref_n = {0.2181529734, 0.4144906495, -0.8835195423};
         unsigned short fv[3] = {1, 3, 2};
         vect n = mapping.normalVector(fv, 0., 0.);
         vect e(p[3], p[2]);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na   << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal TET " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
      {
         vect na;
         std::vector<point> fp = {p[0], p[2], p[3]};
         xmapping::xLagrangeMapping fmapping(xmapping::xReferenceElementType::TRI, fp);
         // mapping.normalVector(&f, 0.,  0., 0. , na);
         const vect ref_n = {-0.1836284555, 0.1159758666, 0.9761302109};
         unsigned short fv[3] = {0, 2, 3};
         vect n = mapping.normalVector(fv, 0., 0.);
         vect e(p[3], p[2]);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na  <<  std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal TET " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
      {
         vect na;
         std::vector<point> fp = {p[0], p[3], p[1]};
         xmapping::xLagrangeMapping fmapping(xmapping::xReferenceElementType::TRI, fp);
         const vect ref_n = {0.2972338858, -0.4564663246, 0.8386241778};
         unsigned short fv[3] = {0, 3, 1};
         vect n = mapping.normalVector(fv, 0., 0.);
         // std::cout << "Normal n  " << std::setprecision(10)<<  na   << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << n << std::endl;
         if (norm(ref_n - n) > eps)
         {
            std::cout << "Error evaluating normal TET " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
   }
   // testing hex case.
   {
      std::cout << "Testing for HEX" << std::endl;
      std::vector<point> p{{0.1, -0.1, 0.2}, {2.1, -0.1, 0.1}, {2.2, 2.1, 0.},   {0.1, 2.05, -0.3},
                           {0., 0., 2.1},    {2., 0., 2.12},   {2.08, 2., 2.13}, {0., 2.01, 2.}};
      xmapping::xLagrangeMapping mapping(xmapping::xReferenceElementType::HEX, p);
      double u, v, w, x = 0., y, z;
      if (mapping.inReferenceElement(-2., 1., 0.))
      {
         std::cout << "Error evaluating inReferenceElement for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (!mapping.inReferenceElement(-0.1, -0.2, 0.))
      {
         std::cout << "Error evaluating inReferenceElement for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      if (!mapping.inReferenceElement(0.3, -0.3, 0.1))
      {
         std::cout << "Error evaluating inReferenceElement for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      };
      mapping.COG(u, v, w);
      if (u != 0. || v != 0. || w != 0.)
      {
         std::cout << "Error evaluating COG for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      mapping.eval(u, v, w, x, y, z);
      if (fabs(x - 1.0725) > eps || fabs(y - 0.995) > eps || fabs(z - 1.04375) > eps)
      {
         std::cout << "Error evaluating eval for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      point min, max;
      point min_ref{0.0, -0.1, -0.3}, max_ref{2.2, 2.1, 2.13};
      mapping.boundingBox(min, max);
      if (norm(min - min_ref) > eps || norm(max - max_ref) > eps)
      {
         std::cout << "Error evaluating boudingBox for HEX" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      const double J_ref = 1.118234188;
      double J_res = mapping.detJac(u, v, w);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating detJac for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens jac_ref{{0.9759723296, -0.004373524844, -0.04127042485},
                   {-0.01695194997, 0.9564465672, 0.07974626514},
                   {0.04921263418, -0.007092655625, 0.9554349276}};
      tens jac_res;
      J_res = mapping.jacInverse(u, v, w, jac_res);
      if (fabs(J_res - J_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for HEX" << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      if (norm(jac_res - jac_ref) > eps)
      {
         std::cout << "Error evaluating jacInverse for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      vect vec_to_push{0.5, 1., 0.3};
      vect vec_pushed_res = vec_to_push;
      vect vec_pushed_ref{0.460486975, 0.9628792947, 0.3146730115};
      mapping.pushBack(0.3, 0.3, 0., 1, &vec_pushed_res);
      if (norm(vec_pushed_res - vec_pushed_ref) > eps)
      {
         std::cout << "Error evaluating PushBack for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      tens tens_to_push{{0.5, 1., 0.}, {3., 2., -0.3}, {-1., 15., 7.}};
      tens tens_pushed_ref{{0.5161360151, 0.3481689072, -0.2875809165},
                           {2.781117462, 3.092135162, 0.2712898858},
                           {-0.9521065774, 14.36655124, 6.69017229}};
      tens tens_pushed_res = tens_to_push;
      mapping.pushBack(0., 0., 0., 1, &tens_pushed_res);
      if (norm(tens_pushed_res - tens_pushed_ref) > eps)
      {
         std::cout << "Error evaluating PushBack for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      double uref = 0.2;
      double vref = -0.3;
      double wref = 0.;
      u = uref;
      v = vref;
      w = wref;
      mapping.eval(u, v, w, x, y, z);
      mapping.invert(x, y, z, u, v, w);
      if (fabs(uref - u) > eps || fabs(vref - v) > eps)
      {
         std::cout << "Error evaluating invert for HEX " << __FILE__ << __LINE__ << std::endl;
         failed = true;
      }
      {
         std::vector<point> fp = {p[0], p[1], p[2], p[3]};
         double u = -0.2;
         double v = 0.1;
         // mapping.normalVector(&f, u,  v, -1. , na);
         const vect ref_n = {0.05965983094, -0.1555910698, -0.9860183181};
         xmapping::xLagrangeMapping fmapping(xmapping::xReferenceElementType::QUAD, fp);
         vect na = -fmapping.normalVector(u, v, 1.);
         unsigned short fv[4] = {0, 1, 2, 3};
         vect nb = mapping.normalVector(fv, u, v);
         unsigned short ifa = 0;
         vect nc = mapping.normalVector(ifa, u, v, -1.);
         std::cout << "Normal n  " << std::setprecision(10) << na << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << nb << std::endl;
         std::cout << "Normal n  " << std::setprecision(10) << nc << std::endl;

         if ((norm(ref_n - na) > eps) || (norm(ref_n - nb) > eps) || (norm(ref_n - nc) > eps))
         {
            std::cout << "Error evaluating normal Hex " << __FILE__ << __LINE__ << std::endl;
            failed = true;
         }
      }
   }
   if (failed) std::cout << "Test Failed " << std::endl;
   return failed;
}
