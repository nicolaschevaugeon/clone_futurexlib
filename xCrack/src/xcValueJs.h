/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _value_sifs__H
#define _value_sifs__H

#include <iostream>
#include <iterator>

#include "mVertex.h"
#include "xValKey.h"
#include "xValue.h"
using namespace std;
using namespace xfem;

namespace xcrack
{
class xcValueJs : public xValue<double>
{
  public:
   xcValueJs() : xValue<double>() { vals.resize(nb_modes, 0.0); }
   double getVal() const override { return vals[mode_ - 1]; }
   void setVal(double in) override
   {
      const bool debug = false;
      if (debug) cout << "mode in setval is " << mode_ << endl;
      vals[mode_ - 1] = in;
      if (debug) cout << "check if value was set " << vals[mode_] << endl;
   }

   static void mode(int i) { mode_ = i; }

   static int getNbModes() { return nb_modes; }
   static void setNbModes(int m) { nb_modes = m; }

   std::ostream& printVal(std::ostream& o) const override
   {
      std::copy(vals.begin(), vals.end(), std::ostream_iterator<double>(o, " "));
      return o;
   }

  private:
   std::vector<double> vals;
   static int nb_modes;
   static int mode_;
};

class xcValueSifs : public xValue<double>
{
  public:
   typedef enum
   {
      GEOM_PLANE_STRAIN,
      GEOM_PLANE_STRESS,
      GEOM_3D
   } geom_t;

   xcValueSifs() : xValue<double>() { vals.resize(nb_modes, 0.0); }
   double getVal() const override { return vals[mode_ - 1]; }
   void setVal(double in) override
   {
      const bool debug = false;
      if (debug) cout << "mode in setval is " << mode_ << endl;
      vals[mode_ - 1] = in;
      if (debug) cout << "check if value was set " << vals[mode_] << endl;
   }

   static void mode(int i) { mode_ = i; }
   static int getNbModes() { return nb_modes; }
   static void setNbModes(int m) { nb_modes = m; }
   static void setYoungAndPoisson(const double& e, const double& p)
   {
      young = e;
      poisson = p;
   }
   static void setGeom(const geom_t& g) { geom = g; }

   std::ostream& printVal(std::ostream& o) const override
   {
      std::copy(vals.begin(), vals.end(), std::ostream_iterator<double>(o, " "));
      return o;
   }

   double getJ1() const
   {
      switch (geom)
      {
         case (GEOM_PLANE_STRESS):
            return (1. / young) * (vals[0] * vals[0] + vals[1] * vals[1]);
            break;
         case (GEOM_PLANE_STRAIN):
            return ((1. - poisson * poisson) / young) * (vals[0] * vals[0] + vals[1] * vals[1]);
            break;
         case (GEOM_3D):
            return ((1. - poisson * poisson) / young) * (vals[0] * vals[0] + vals[1] * vals[1]) +
                   ((1. + poisson) / young) * vals[2] * vals[2];
            break;
         default:
            throw;
            return 0.;
      }
   }

  private:
   std::vector<double> vals;
   static int nb_modes;
   static int mode_;
   static geom_t geom;
   static double young;
   static double poisson;
};

class xcValueJsAndSifs : public xValue<double>
{
   typedef enum
   {
      JS,
      SIFS,
      ERR
   } field_t;

  public:
   double getVal() const override
   {
      switch (field)
      {
         case (JS):
            return Js.getVal();
            break;
         case (SIFS):
            return Sifs.getVal();
            break;
         case (ERR):
            return compute_err();
            break;
         default:
            throw;
            return 0.;
      }
   }
   void setVal(double in) override
   {
      switch (field)
      {
         case (JS):
            Js.setVal(in);
            break;
         case (SIFS):
            Sifs.setVal(in);
            break;
         case (ERR):
            break;
         default:
            throw;
      }
   }

   static void js() { field = JS; }
   static void sifs() { field = SIFS; }
   static void js(int i)
   {
      field = JS;
      xcValueJs::mode(i);
   }
   static void sifs(int i)
   {
      field = SIFS;
      xcValueSifs::mode(i);
   }
   static void err() { field = ERR; }

   std::ostream& print(std::ostream& o) const override
   {
      o << "Js  " << endl;
      Js.print(o) << endl;
      o << "Sifs " << endl;
      Sifs.print(o) << endl;
      return o;
   }

   std::ostream& printVal(std::ostream& o) const override { return o; };

  protected:
   double compute_err() const
   {
      js(1);
      double j = getVal();
      double j_from_K = Sifs.getJ1();
      double errK = (j > 0.) ? (sqrt(j) - sqrt(j_from_K)) : (-sqrt(-j) - sqrt(j_from_K));
      if (j_from_K > 1.e-16) errK /= sqrt(j_from_K);
      err();
      return errK;
   }

  private:
   static field_t field;
   xcValueJs Js;
   xcValueSifs Sifs;
};

}  // end namespace xcrack
#endif
