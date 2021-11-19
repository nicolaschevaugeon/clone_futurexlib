/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _ENV_H
#define _ENV_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <string>

#include "xDebug.h"
#include "xEval.h"
#include "xGeomElem.h"
#include "xStringManager.h"
#include "xTensor2.h"
#include "xVector.h"

namespace xfem
{
// possible value of the Geom member of xEnv
#define BC_POINT 0
#define BC_LINE 1
#define BC_SURFACE 2
#define BC_VOLUME 3
#define ZONE_ENV 4
#define GROUP_INT 5
#define BC_POINT_NAME 6
#define BC_LINE_NAME 7
#define BC_SURFACE_NAME 8
#define BC_VOLUME_NAME 9
#define IC_POINT 10
#define IC_LINE 11
#define IC_SURFACE 12
#define IC_VOLUME 13
#define RIGID_LINE 14
#define RIGID_SURFACE 15
#define CONTACT_LINE 16

// The order is very important because the priority if fix
// possible value of the Type member of xEnv
#define FIX 1
#define ASS 2
#define FIX_AND_MEASURE 3

class xPieceWiseLinearDoubleToDouble;

typedef xtool::xSingleton<xtool::xStringManager<int>> xNamesSingleton;

class xEnv
{
  public:
   std::string Phys;       // key1
   int Geom;               // key2
   int Entity;             // key3
   std::string Geom_name;  // key4
   int Type;
   double Val_fix;
   int Val_ass;

  private:
   xPieceWiseLinearDoubleToDouble* evolution;

  public:
   xtool::xStringManager<int>& names;

   int getDimension() const;

   ~xEnv();
   xEnv();
   // PHYS GEOM ENTITY
   xEnv(const std::string& key1, const int& key2, const int& key3);
   xEnv(const std::string& key1, const int& key2, const std::string& key4);
   xEnv(const xEnv& rhs);
   void clear();

   friend bool operator<(const xEnv& c1, const xEnv& c2) { return (compare(c1, c2) == -1) ? true : false; }
   friend bool operator>(const xEnv& c1, const xEnv& c2) { return (compare(c1, c2) == 1) ? true : false; }
   friend bool operator==(const xEnv& c1, const xEnv& c2) { return (compare(c1, c2) == 0) ? true : false; }
   friend bool operator!=(const xEnv& c1, const xEnv& c2) { return (compare(c1, c2) != 0) ? true : false; }

   static int compare(const xEnv& c1, const xEnv& c2);
   void print(std::ostream& fp) const;
   void print() const;

   void defineFixed(const std::string& key1, int key2, int key3, double val);
   void defineFixedAndMeas(const std::string& key1, int key2, int key3, double val);
   void defineAssociate(const std::string& key1, int key2, int key3, int ass);

   double getValue() const { return Val_fix; }

   void setEvolution(const std::map<double, double>& evo);

   const xPieceWiseLinearDoubleToDouble& getEvolution() const;

   // void setHomogeneous(void) {Val_fix = 0.0;}

   bool hasAnEvolution() const { return (evolution != nullptr); }

   // celine begin
   std::string getGeomName() const { return Geom_name; }
   void set(const int& key2, const int& key3);
   void set(const std::string& key4, const int& key2);
   void defineFixed(const std::string& key1, double val);
   void defineAssociate(const std::string& key1, int val);
   // celine end
};

class xPieceWiseLinearDoubleToDouble
{
  public:
   xPieceWiseLinearDoubleToDouble(const std::map<double, double>& evo);
   double operator()(const double& t) const;

  private:
   std::map<double, double> evolution;
   typedef std::map<double, double>::const_iterator const_iterator;
   typedef std::map<double, double>::value_type value_type;
   double tmin, tmax;
};

class xEvalNormal : public xEval<xtensor::xVector<>>
{
  public:
   xEvalNormal();
   xEvalNormal(const xfem::xEntityToEntity& upper_);
   void operator()(const xGeomElem* appro, const xGeomElem* integ, xtensor::xVector<>&) const override;

  private:
   xfem::xEntityToEntity upper;
};

class xEvalLinearInTimeVector : public xEval<xtensor::xVector<>>
{
  public:
   xEvalLinearInTimeVector(const xtensor::xVector<>& vfin, double tf) : val_final(vfin), t_final(tf) { assert(t_final > 0.); }
   void setTime(const double& t_) { t = t_; }
   void operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& result) const override
   {
      const bool debug = xdebug_flag;
      if (debug) std::cout << " t_final " << t_final << " t " << t << std::endl;
      if (debug) std::cout << "v_final " << val_final << std::endl;
      if (debug) std::cout << "result " << val_final * (t / t_final) << std::endl;
      result = val_final * (t / t_final);
   }

  private:
   xEvalLinearInTimeVector();
   xtensor::xVector<> val_final;
   double t_final, t;
};

class xEvalVectFromPoint : public xEval<xtensor::xVector<>>
{
  public:
   xEvalVectFromPoint(const xtensor::xPoint& orig_);
   void operator()(const xGeomElem* appro, const xGeomElem* integ, xtensor::xVector<>&) const override;

  private:
   const xtensor::xPoint orig;
};

class xPhysicalEnv
{
  public:
   typedef std::set<xEnv> stl_container_t;
   typedef stl_container_t::const_iterator const_iterator;
   typedef stl_container_t::iterator iterator;

   iterator begin() { return info.begin(); }
   iterator end() { return info.end(); }
   const_iterator begin() const { return info.begin(); }
   const_iterator end() const { return info.end(); }
   int size() const { return info.size(); }
   void clear() { info.clear(); }
   const_iterator find(const xEnv& d) const { return info.find(d); }

   void print(std::ostream& fp)
   {
      fp << "The total number of env info is " << info.size() << std::endl;
      for (const_iterator it = info.begin(); it != info.end(); it++)
      {
         it->print(fp);
      }
      return;
   }

   void add(const xEnv& env)
   {
      if (find(env) != end())
      {
         printf("Info: An already existing element is being added\n");
         printf("Check the < operator for Info_c to understand why\n");
         printf("or verify if a Boundary Condition is not multiple defined\n");
         // assert(1 == 0);
         // assert removed because a ligne or a point may
         // be associted to more than one point
      }
      info.insert(env);
   }

  private:
   stl_container_t info;
};

// class xPhysicalEnv : public xPhysicalEnvBase{

// public:
//    xPhysicalEnv(void) : xPhysicalEnvBase(), HomogeneousModeOn(false) {}
//    void SetHomogeneousModeOn(void)  {HomogeneousModeOn = true;}
//    void SetHomogeneousModeOff(void) {HomogeneousModeOn = false;}
//    virtual xEnv currentItem() const {
//      xEnv env = *itcurr;
//      if (HomogeneousModeOn) {env.setHomogeneous();}
//      return env;
//    }

// private:
//  bool HomogeneousModeOn;

//  };

}  // namespace xfem

#endif
