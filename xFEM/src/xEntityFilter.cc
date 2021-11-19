/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#include "xEntityFilter.h"

#include <iostream>
#include <string>

#include "GEntity.h"
#include "mEntity.h"
#include "mVertex.h"
#include "xDebug.h"
#include "xDomain.h"
#include "xLevelSet.h"
#include "xMesh.h"
#include "xSubMesh.h"
#include "xZone.h"

using namespace std;
using AOMD::mEntity;
using AOMD::mVertex;

namespace xfem
{
xAcceptSupport::xAcceptSupport(xLevelSet& ls, xMesh* m) : LS(&ls), mesh(m) {}

bool xAcceptSupport::operator()(AOMD::mEntity* e) const
{
   std::set<mEntity*> support_e;
   mesh->lookupSupport(e, support_e);

   std::set<mEntity*>::iterator itc;
   std::vector<double> vals;
   bool test = false;

   for (itc = support_e.begin(); itc != support_e.end(); ++itc)
   {
      AOMD::mEntity* es = *itc;
      vals.clear();
      vals = LS->getVals(es);
      if (*min_element(vals.begin(), vals.end()) <= 0.0)
      {
         test = true;
         break;
      }
   }
   return test;
}

xAcceptCutLevelSet::xAcceptCutLevelSet(xLevelSet& ls, double _dist, xEntityToEntity _ent) : LS(&ls), dist(_dist), ent(_ent) {}

bool xAcceptCutLevelSet::operator()(AOMD::mEntity* ei) const
{
   mEntity* e = ent(ei);
   std::vector<double> vals = LS->getVals(e);
   bool test = false;
   double ls_max = *max_element(vals.begin(), vals.end());
   double ls_min = *min_element(vals.begin(), vals.end());

   if (ls_max * ls_min <= dist) test = true;
   return test;
}

xAcceptStrictCutLevelSet::xAcceptStrictCutLevelSet(xLevelSet& ls, double _dist, xEntityToEntity _ent)
    : LS(&ls), dist(_dist), ent(_ent)
{
}

bool xAcceptStrictCutLevelSet::operator()(AOMD::mEntity* ei) const
{
   mEntity* e = ent(ei);
   std::vector<double> vals = LS->getVals(e);
   bool test = false;
   double ls_max = *max_element(vals.begin(), vals.end());
   double ls_min = *min_element(vals.begin(), vals.end());

   if (ls_max * ls_min < dist) test = true;
   return test;
}

xAcceptStrictlyInsideLevelSet::xAcceptStrictlyInsideLevelSet(xLevelSet& ls, double _dist, xEntityToEntity _ent)
    : LS(&ls), dist(_dist), ent(_ent)
{
}

bool xAcceptStrictlyInsideLevelSet::operator()(AOMD::mEntity* ei) const
{
   mEntity* e = ent(ei);
   std::vector<double> vals = LS->getVals(e);
   bool test = false;
   if (*max_element(vals.begin(), vals.end()) <= dist) test = true;
   return test;
}

xAcceptInsideAndCrossingLevelSet::xAcceptInsideAndCrossingLevelSet(xLevelSet& ls, xEntityToEntity _ent) : LS(&ls), ent(_ent) {}

bool xAcceptInsideAndCrossingLevelSet::operator()(AOMD::mEntity* ei) const
{
   mEntity* e = ent(ei);
   std::vector<double> vals = LS->getVals(e);
   bool test = false;
   if (*min_element(vals.begin(), vals.end()) <= 0.0) test = true;
   return test;
}

xAcceptCurvatureLessThan::xAcceptCurvatureLessThan(xLevelSet& ls, double _curv, xEntityToEntity _ent)
    : LS(&ls), curv(_curv), ent(_ent)
{
}

bool xAcceptCurvatureLessThan::operator()(AOMD::mEntity* ei) const
{
   mEntity* e = ent(ei);
   double vals = (LS->getTrueCurv(e)).trace();
   return (fabs(1. / vals) >= curv);
}

xAccept::xAccept(const std::string& z) : zone_id{xZone::names.getId(z)} {}

xAccept::xAccept(const int& zone_id_) : zone_id{zone_id_} {}

bool xAccept::operator()(AOMD::mEntity* e) const
{
   const bool debug = false;
   int* pid = xDomain::get(*e);
   int id = (pid) ? (*pid) : e->getClassification()->tag();
   if (debug) cout << "in EntityFilter accept id is " << id << endl;
   return id == zone_id;
}

xAcceptUnion::xAcceptUnion(const std::list<xEntityFilter>& _filters) : filtersList(_filters) {}
bool xAcceptUnion::operator()(AOMD::mEntity* e) const
{
   std::list<xEntityFilter>::const_iterator itb = filtersList.begin();
   std::list<xEntityFilter>::const_iterator ite = filtersList.end();
   while (itb != ite)
   {
      if ((*itb)(e)) return true;
      ++itb;
   }
   return false;
}

xAcceptUnionBinary::xAcceptUnionBinary(const xEntityFilter f1_, const xEntityFilter f2_) : f1(f1_), f2(f2_) {}
bool xAcceptUnionBinary::operator()(AOMD::mEntity* e) const { return (f1(e) || f2(e)); }

xAcceptIntersection::xAcceptIntersection(const std::list<xEntityFilter>& _filters) : filtersList(_filters) {}
bool xAcceptIntersection::operator()(AOMD::mEntity* e) const
{
   const bool debug = xdebug_flag;
   std::list<xEntityFilter>::const_iterator itb = filtersList.begin();
   std::list<xEntityFilter>::const_iterator ite = filtersList.end();
   while (itb != ite)
   {
      if (!((*itb)(e))) return false;
      ++itb;
   }
   if (debug) std::cout << "true !!!" << std::endl;
   return true;
}

xAcceptIntersectionBinary::xAcceptIntersectionBinary(const xEntityFilter f1_, const xEntityFilter f2_) : f1(f1_), f2(f2_) {}
bool xAcceptIntersectionBinary::operator()(AOMD::mEntity* e) const { return (f1(e) && f2(e)); }

xAcceptIntersectionTernary::xAcceptIntersectionTernary(const xEntityFilter f1_, const xEntityFilter f2_, const xEntityFilter f3_)
    : f1(f1_), f2(f2_), f3(f3_)
{
}
bool xAcceptIntersectionTernary::operator()(AOMD::mEntity* e) const { return (f1(e) && f2(e) && f3(e)); }

bool xAcceptUpperAdjacency::operator()(AOMD::mEntity* e) const
{
   const bool debug = xdebug_flag;
   if (debug)
   {
      cout << " in xAcceptUpperAdjacency " << endl;
      cout << " e " << endl;
      e->print();
   }
   mEntity* e_upper = e->get(e->getLevel() + 1, 0);
   if (debug)
   {
      assert(e_upper != nullptr);
      cout << " e_upper " << endl;
      e_upper->print();
   }
   return accept(e_upper);
}

xAcceptOnBoundaryOfRegion::xAcceptOnBoundaryOfRegion(const xRegion& reg_) : region(reg_) {}
bool xAcceptOnBoundaryOfRegion::operator()(AOMD::mEntity* e) const
{
   assert(e->getLevel() < 3);
   mEntity* upper1 = e->get(e->getLevel() + 1, 0);
   bool test1 = region.IsInRegion(upper1);
   bool test2 = false;
   if (e->size(e->getLevel() + 1) == 2)
   {
      mEntity* upper2 = e->get(e->getLevel() + 1, 1);
      test2 = region.IsInRegion(upper2);
   }
   if (test1 && test2) return false;
   if (!test1 && !test2) return false;
   return true;
}

xAcceptOnBoundaryOfSubMesh::xAcceptOnBoundaryOfSubMesh(const xSubMesh& _sub) : sub(_sub) {}
bool xAcceptOnBoundaryOfSubMesh::operator()(AOMD::mEntity* e) const { return sub.isOnBoundary(e); }

xAcceptInMesh::xAcceptInMesh(const xMesh& _mesh) : mesh(_mesh) {}

bool xAcceptInMesh::operator()(AOMD::mEntity* e) const { return mesh.find(e); }

xAcceptInSubMesh::xAcceptInSubMesh(const xSubMesh& _sub) : sub(_sub) {}

bool xAcceptInSubMesh::operator()(AOMD::mEntity* e) const { return sub.find(e); }

xAcceptLessThanForOneNode::xAcceptLessThanForOneNode(std::function<double(AOMD::mVertex*)>& _eval, double _limit)
    : eval(_eval), limit(_limit)
{
}

bool xAcceptLessThanForOneNode::operator()(AOMD::mEntity* e) const
{
   int size0 = e->size(0);
   for (int k = 0; k < size0; ++k)
   {
      double val = eval(dynamic_cast<AOMD::mVertex*>(e->get(0, k)));
      if (val <= limit) return true;
   }
   return false;
}

bool xAcceptOnBoundary::operator()(AOMD::mEntity* e) const
{
   assert(e->getLevel() < 3);
   if (e->size(e->getLevel() + 1) == 1) return true;
   return false;
}

bool xAcceptOnClassifiedBoundary::operator()(AOMD::mEntity* e) const
{
   int level = e->getLevel();
   assert(level < 3);
   if (GEN_type(e->getClassification()) == level) return true;
   return false;
}
xAcceptOnClassification::xAcceptOnClassification(const int iclass_, int iwhat_) : iclass(iclass_), iwhat(iwhat_) {}
bool xAcceptOnClassification::operator()(AOMD::mEntity* e) const
{
   pGEntity g = e->getClassification();
   if (!e->isAdjacencyCreated(e->getLevel()))
      if (GEN_tag(g) == iclass && GEN_type(g) == iwhat) return true;
   return false;
}

xAcceptInVolume::xAcceptInVolume(std::function<bool(const double& x, const double& y, const double& z)>& in_vol_)
    : in_vol(in_vol_)
{
}

bool xAcceptInVolume::operator()(AOMD::mEntity* e) const
{
   if (e->getLevel())
   {
      for (int i = 0; i < e->size(0); ++i)
      {
         xtensor::xPoint p(((mVertex*)e->get(0, i))->point());
         if (!in_vol(p(0), p(1), p(2))) return false;
      }
      return true;
   }
   else
   {
      xtensor::xPoint p(((mVertex*)e)->point());
      return in_vol(p(0), p(1), p(2));
   }
}

xAcceptOnIntTag::xAcceptOnIntTag(const unsigned int tag_, const int accepted_tag_val_)
    : tag(tag_), accepted_tag_val(accepted_tag_val_)
{
}

bool xAcceptOnIntTag::operator()(AOMD::mEntity* e) const
{
   if (e->getAttachedInt(tag) == accepted_tag_val) return true;
   return false;
}

}  // namespace xfem
