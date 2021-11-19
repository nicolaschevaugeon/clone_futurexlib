#ifndef _XPOINT_H_
#define _XPOINT_H_

#include <math.h>

#include <iostream>

namespace xtensor
{
class xPoint
{
   double pos[3];

  public:
   inline xPoint(double x = 0.0, double y = 0.0, double z = 0.0)
   {
      pos[0] = x;
      pos[1] = y;
      pos[2] = z;
   }

   template <typename T>
   inline xPoint(T other)
   {
      pos[0] = other(0);
      pos[1] = other(1);
      pos[2] = other(2);
   }

   inline bool lexicographicLessThan(const xPoint& other, double EPS) const
   {
      if (pos[0] < other(0) - EPS) return 1;
      if (pos[0] > other(0) + EPS) return 0;
      if (pos[1] < other(1) - EPS) return 1;
      if (pos[1] > other(1) + EPS) return 0;
      if (pos[2] < other(2) - EPS) return 1;
      if (pos[2] > other(2) + EPS) return 0;
      return 0;
   }

   friend std::ostream& operator<<(std::ostream& out, const xPoint& obj)
   {
      out << "(" << obj(0) << "," << obj(1) << "," << obj(2) << ")";
      return out;
   }

   void print() { std::cout << " xPoint=" << pos[0] << "," << pos[1] << "," << pos[2] << std::endl; }

   inline double& operator()(int i)
   {
      if (i >= 3) throw -1;
      return pos[i];
   }

   inline double operator()(int i) const
   {
      if (i >= 3) throw -1;
      return pos[i];
   }

   inline double& operator[](int i)
   {
      if (i >= 3) throw -1;
      return pos[i];
   }

   inline double operator[](int i) const
   {
      if (i >= 3) throw -1;
      return pos[i];
   }

   inline xPoint operator*(double other) const { return xPoint(pos[0] * other, pos[1] * other, pos[2] * other); }
   inline xPoint operator/(double other) const { return xPoint(pos[0] / other, pos[1] / other, pos[2] / other); }

   inline xPoint& operator+=(const xPoint& other)
   {
      pos[0] += other.pos[0];
      pos[1] += other.pos[1];
      pos[2] += other.pos[2];
      return *this;
   }
   inline xPoint operator+(const xPoint& other) const
   {
      return xPoint(pos[0] + other.pos[0], pos[1] + other.pos[1], pos[2] + other.pos[2]);
   }

   inline xPoint& operator*=(double other)
   {
      pos[0] *= other;
      pos[1] *= other;
      pos[2] *= other;
      return *this;
   }
   inline xPoint& operator/=(double other)
   {
      pos[0] /= other;
      pos[1] /= other;
      pos[2] /= other;
      return *this;
   }
   inline xPoint& operator-=(const xPoint& other)
   {
      pos[0] -= other.pos[0];
      pos[1] -= other.pos[1];
      pos[2] -= other.pos[2];
      return *this;
   }

   inline xPoint operator-(const xPoint& other) const
   {
      return xPoint(pos[0] - other.pos[0], pos[1] - other.pos[1], pos[2] - other.pos[2]);
   }

   inline bool operator<(const xPoint& other) const { return lexicographicLessThan(other, 1.e-6); }

   inline bool closeTo(const xPoint& other, double EPS) const
   {
      return ((pow(pos[0] - other(0), 2) + pow(pos[1] - other(1), 2) + pow(pos[2] < other(2), 2)) < pow(EPS, 2));
   }
   inline bool operator==(const xPoint& other) const
   {
      return (pos[0] == other.pos[0] && pos[1] == other.pos[1] && pos[2] == other.pos[2]);
   }

};  // end of xPoint

}  // namespace xtensor

#endif
