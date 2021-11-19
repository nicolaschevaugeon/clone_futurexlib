/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#include "xEnv.h"

#include <cassert>
#include <iostream>
#include <string>

#include "mAOMD.h"
#include "mEntity.h"
#include "xDebug.h"
#include "xGeomElem.h"
#include "xMesh.h"
#include "xStringManager.h"
#include "xTensor2.h"
#include "xVector.h"

namespace xfem
{
using AOMD::mEntity;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

xEnv::~xEnv()
{
   if (evolution) delete evolution;
}

// Constructors
xEnv::xEnv()
    : Phys(""),
      Geom(0),
      Entity(0),
      Geom_name(""),
      Type(0),
      Val_fix(0.0),
      Val_ass(0),
      evolution(nullptr),
      names(xNamesSingleton::instance())
{
}
xEnv::xEnv(const string& key1, const int& key2, const int& key3)
    : Phys(key1),
      Geom(key2),
      Entity(key3),
      Geom_name(xNamesSingleton::instance().getName(key3)),
      Type(0),
      Val_fix(0.0),
      Val_ass(0),
      evolution(nullptr),
      names(xNamesSingleton::instance())
{
}

xEnv::xEnv(const string& key1, const int& key2, const string& key4)
    : Phys(key1),
      Geom(key2),
      Entity(xNamesSingleton::instance().getId(key4)),
      Geom_name(key4),
      Type(0),
      Val_fix(0.0),
      Val_ass(0),
      evolution(nullptr),
      names(xNamesSingleton::instance())
{
}

void xEnv::set(const int& key2, const int& key3)
{
   Geom = key2;
   Entity = key3;
   Geom_name = names.getName(key3);
   Type = 0;
   Val_fix = 0.0;
   Val_ass = 0;
   evolution = nullptr;
}

void xEnv::set(const string& key4, const int& key2)
{
   Geom = key2;
   Entity = names.getId(key4);
   Geom_name = key4;
   Type = 0;
   Val_fix = 0.0;
   Val_ass = 0;
   evolution = nullptr;
}

xEnv::xEnv(const xEnv& rhs) : names(rhs.names)
{
   Phys = rhs.Phys;
   Geom = rhs.Geom;
   Entity = rhs.Entity;
   Geom_name = rhs.Geom_name;
   Type = rhs.Type;
   Val_fix = rhs.Val_fix;
   Val_ass = rhs.Val_ass;
   if (rhs.evolution)
      evolution = new xPieceWiseLinearDoubleToDouble(*rhs.evolution);
   else
      evolution = nullptr;
}

void xEnv::clear()
{
   delete evolution;
   evolution = nullptr;
}

int xEnv::compare(const xEnv& c1, const xEnv& c2)
{
   if (c1.Entity > c2.Entity) return 1;
   if (c1.Entity < c2.Entity) return -1;
   if (c1.Phys > c2.Phys) return 1;
   if (c1.Phys < c2.Phys) return -1;
   if (c1.Geom > c2.Geom) return 1;
   if (c1.Geom < c2.Geom) return -1;
   if (c1.Geom_name > c2.Geom_name) return 1;
   if (c1.Geom_name < c2.Geom_name) return -1;
   return 0;
}

void xEnv ::print(std::ostream& fp) const
{
   if (Type == FIX)
   {
      fp << Val_fix << endl;
      return;
   }
   else if (Type == ASS)
   {
      fp << Val_ass << endl;
      return;
   }
   else if (Type == FIX_AND_MEASURE)
   {
      fp << Val_fix << endl;
      return;
   }
   return;
}

void xEnv ::print() const
{
   print(cout);
   return;
}

void xEnv ::defineFixed(const string& key1, int key2, int key3, double val)
{
   Phys = key1;
   Geom = key2;
   Entity = key3;
   Type = FIX;
   Val_fix = val;
}

void xEnv ::defineFixed(const string& key1, double val)
{
   Phys = key1;
   Type = FIX;
   Val_fix = val;
}

void xEnv ::defineFixedAndMeas(const string& key1, int key2, int key3, double val)
{
   Phys = key1;
   Geom = key2;
   Entity = key3;
   Type = FIX_AND_MEASURE;
   Val_fix = val;
}

void xEnv ::defineAssociate(const string& key1, int key2, int key3, int val)
{
   Phys = key1;
   Geom = key2;
   Entity = key3;
   Type = ASS;
   Val_ass = val;
}

void xEnv ::defineAssociate(const string& key1, int val)
{
   Phys = key1;
   Type = ASS;
   Val_ass = val;
}

int xEnv::getDimension() const
{
   int dim;
   switch (Geom)
   {
      case BC_POINT:
         dim = 0;
         break;
      case BC_POINT_NAME:
         dim = 0;
         break;
      case BC_LINE:
         dim = 1;
         break;
      case BC_LINE_NAME:
         dim = 1;
         break;
      case BC_SURFACE:
         dim = 2;
         break;
      case BC_SURFACE_NAME:
         dim = 2;
         break;
      case BC_VOLUME:
         dim = 3;
         break;
      case IC_POINT:
         dim = 0;
         break;
      case IC_LINE:
         dim = 1;
         break;
      case IC_SURFACE:
         dim = 2;
         break;
      case IC_VOLUME:
         dim = 3;
         break;
      case RIGID_LINE:
         dim = 1;
         break;
      case RIGID_SURFACE:
         dim = 2;
         break;
      case CONTACT_LINE:
         dim = 1;
         break;
      default:
         assert(0);
         throw;
         break;
   }
   return dim;
}

void xEnv::setEvolution(const std::map<double, double>& evo) { evolution = new xPieceWiseLinearDoubleToDouble(evo); }

const xPieceWiseLinearDoubleToDouble& xEnv::getEvolution() const { return *evolution; }

xPieceWiseLinearDoubleToDouble::xPieceWiseLinearDoubleToDouble(const std::map<double, double>& evo) : evolution(evo)
{
   const bool debug = xdebug_flag;
   if (evolution.size() <= 1)
   {
      std::cerr << " wrong evolution, only " << evolution.size() << " points " << std::endl;
      assert(0);
   }
   tmin = evo.begin()->first;
   tmax = (--evo.end())->first;
   if (debug)
   {
      cout << "evolution" << endl;
      for (const_iterator it = evolution.begin(); it != evolution.end(); ++it) cout << it->first << " " << it->second << endl;
   }
   if (debug) cout << " In xdynPieceWiseLinearDoubleToDouble, tmin is " << tmin << " tmax is " << tmax << endl;
}

double xPieceWiseLinearDoubleToDouble::operator()(const double& t) const
{
   const bool debug = xdebug_flag;
   if (t < tmin || t > tmax)
   {
      cerr << " time t = " << t << " is outside the bounds [ " << tmin << " , " << tmax << "]" << endl;
      assert(0);
   }
   if (debug) cout << "looking for PieceWiseLinear function  for value " << t << endl;
   double val;
   const_iterator lower = evolution.lower_bound(t);
   const_iterator upper = evolution.upper_bound(t);
   if (lower == upper) lower--;
   val = lower->second + (t - lower->first) / (upper->first - lower->first) * (upper->second - lower->second);
   if (debug) cout << " t is " << t << " and val is " << val << endl;
   return val;
}

xEvalNormal::xEvalNormal()
    : upper([](AOMD::mEntity* e) {
         mEntity* const* pup = xMesh::get_const_is_in_partition_of().getData(*e);
         return pup ? *pup : nullptr;
         // return e->getAttachedEntity(xMesh::get_is_in_partition_of_tag());
      })
{
}
xEvalNormal::xEvalNormal(const xfem::xEntityToEntity& upper_) : upper(upper_) {}

void xEvalNormal::operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& normal) const
{
   const bool debug = xdebug_flag;
   if (debug)
   {
      cout << " bef normal computation " << endl;
      xtensor::xPoint p = geo_integ->getXYZ();
      cout << " current point " << p << endl;
      cout << "appro " << endl;
      geo_appro->getEntity()->print();
      cout << "integ " << endl;
      geo_integ->getEntity()->print();
   }
   mEntity* cur = geo_integ->getEntity();
   mEntity* up = upper(cur);
   if (up)
   {
      if (up->getLevel() == cur->getLevel())
      {
         geo_appro->normalVector(up, normal);
      }
      else
      {
         geo_appro->normalVector(cur, normal);
      }
   }
   else
   {
      geo_appro->normalVector(cur, normal);
   }
   if (debug) cout << " normal is  " << normal << endl;
}

xEvalVectFromPoint::xEvalVectFromPoint(const xtensor::xPoint& orig_) : orig(orig_) {}

void xEvalVectFromPoint::operator()(const xGeomElem* geo_appro, const xGeomElem* geo_integ, xtensor::xVector<>& vec) const
{
   vec = xtensor::xVector<>(orig, geo_integ->getXYZ());
}

}  // namespace xfem
