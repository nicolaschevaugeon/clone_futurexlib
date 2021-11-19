/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.

*/

#ifndef _xcFilters_
#define _xcFilters_
#include "xAOMDEntityUtil.h"
#include "xLevelSet.h"
#include "xPointToDouble.h"

struct xcAcceptCylinderAroundTipLevelSet
{
   xcAcceptCylinderAroundTipLevelSet(const xfem::xLevelSet& _lsn, const xfem::xLevelSet& _lst, double _radius)
       : lsn(_lsn), lst(_lst), radius(_radius)
   {
   }
   bool operator()(AOMD::mEntity* e)
   {
      auto nds = xinterface::aomd::getVertices(*e);
      for (auto v : nds)
      {
         double X = lsn(v);
         double Y = lst(v);
         double rr = sqrt(X * X + Y * Y);
         if (rr < radius) return true;
      }
      return false;
   }

  private:
   const xfem::xLevelSet& lsn;
   const xfem::xLevelSet& lst;
   double radius;
};

// return true if at least one node is inside radius.
struct xcAcceptCylinderAroundTip
{
   xcAcceptCylinderAroundTip(std::pair<xfem::xPointToDouble, xfem::xPointToDouble> _XY, double rad) : XY(_XY), r(rad) {}
   bool operator()(AOMD::mEntity* e)
   {
      auto nds = xinterface::aomd::getPoints(*e);
      for (auto P : nds)
      {
         double X = (XY.first)(P);
         double Y = (XY.second)(P);
         double rr = sqrt(X * X + Y * Y);
         if (rr < r) return true;
      }
      return false;
   }
   std::pair<xfem::xPointToDouble, xfem::xPointToDouble> XY;
   double r;
};

class xcAcceptCyl
{
  public:
   xcAcceptCyl(double rad) { r = rad; }
   bool operator()(AOMD::mEntity* e)
   {
      auto nds = xinterface::aomd::getPoints(*e);
      for (auto P : nds)
      {
         double rr = sqrt(P(0) * P(0) + P(1) * P(1));
         if (rr < r) return true;
      }
      return false;
   }

  private:
   double r;
};

class xAcceptSquare
{
  public:
   xAcceptSquare(double rad) : r{rad} {}
   bool operator()(AOMD::mEntity* e)
   {
      auto nds = xinterface::aomd::getPoints(*e);
      for (auto P : nds)
      {
         double rr = fabs(P(0));
         if (fabs(P(1)) > rr) rr = fabs(P(1));
         if (rr < r) return true;
      }
      return false;
   }

  private:
   double r;
};
#endif
