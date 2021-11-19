/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef __ENTITY_FILTER__H
#define __ENTITY_FILTER__H

#include <functional>
#include <list>
#include <string>
#include "xGenericOperations.h"

namespace AOMD
{
class mEntity;
class mVertex;
}  // namespace AOMD

namespace xfem
{
typedef std::function<AOMD::mEntity*(AOMD::mEntity*)> xEntityToEntity;
class xMesh;
class xLevelSet;
class xRegion;
class xSubMesh;

typedef std::function<bool(AOMD::mEntity* e)> xEntityFilter;

//! entityFilter, that always return true
struct xAcceptAll
{
   bool operator()(AOMD::mEntity* e) const { return true; }
};

//! entityFilter, that always return false
struct xAcceptNone
{
   bool operator()(AOMD::mEntity* e) const { return false; }
};

//! entityFilter given a list of filter return true if one filter is true.
struct xAcceptUnion
{
  public:
   xAcceptUnion(const std::list<xEntityFilter>&);
   bool operator()(AOMD::mEntity* e) const;

  private:
   std::list<xEntityFilter> filtersList;
};

//! entityFilter given 2 filter, return true if one of the filter is true.
//! this is a particular case of xAcceptUnion implemented for performence purpose (no iterator, loop and so on)
class xAcceptUnionBinary
{
  public:
   xAcceptUnionBinary(const xEntityFilter f1_, const xEntityFilter f2_);
   bool operator()(AOMD::mEntity* e) const;

  private:
   xEntityFilter f1;
   xEntityFilter f2;
};

struct xAcceptOnBoundaryOfSubMesh
{
  public:
   xAcceptOnBoundaryOfSubMesh(const xSubMesh& sub);
   bool operator()(AOMD::mEntity* e) const;

  private:
   const xSubMesh& sub;
};

struct xAcceptInMesh
{
  public:
   xAcceptInMesh(const xMesh& sub);
   bool operator()(AOMD::mEntity* e) const;

  private:
   const xMesh& mesh;
};

struct xAcceptInSubMesh
{
  public:
   xAcceptInSubMesh(const xSubMesh& sub);
   bool operator()(AOMD::mEntity* e) const;

  private:
   const xSubMesh& sub;
};

//! entityFilter given a list of filter, return true if all filter are true.
struct xAcceptIntersection
{
  public:
   xAcceptIntersection(const std::list<xEntityFilter>&);
   bool operator()(AOMD::mEntity* e) const;

  private:
   std::list<xEntityFilter> filtersList;
};

//! entityFilter given 2 filter, return true if all filter are true.
//! this is a particular case of xAcceptIntersection implemented for performence purpose (no iterator, loop and so on)
class xAcceptIntersectionBinary
{
  public:
   xAcceptIntersectionBinary(const xEntityFilter f1_, const xEntityFilter f2_);
   bool operator()(AOMD::mEntity* e) const;

  private:
   xEntityFilter f1;
   xEntityFilter f2;
};

//! entityFilter given 3 filters, return true if all filter are true.
//! this is a particular case of xAcceptIntersection implemented for performence purpose (no iterator, loop and so on)
class xAcceptIntersectionTernary
{
  public:
   xAcceptIntersectionTernary(const xEntityFilter f1_, const xEntityFilter f2_, const xEntityFilter f3_);
   bool operator()(AOMD::mEntity* e) const;

  private:
   xEntityFilter f1;
   xEntityFilter f2;
   xEntityFilter f3;
};

//! entityFilter, given an filter at construction, return the contrary of the input filter
struct xAcceptInvert
{
  public:
   xAcceptInvert(const xEntityFilter& f) : filter(f) {}
   bool operator()(AOMD::mEntity* e) const { return !filter(e); }

  private:
   xEntityFilter filter;
};

//! entityFilter, return true if the entity is in the support of the levelset
struct xAcceptSupport
{
  public:
   xAcceptSupport(xLevelSet& ls, xMesh* m);
   bool operator()(AOMD::mEntity* e) const;

  private:
   xLevelSet* LS;
   xMesh* mesh;
};

//! entityFilter, return true if the entity is cut by the iso-zero of the level-set
struct xAcceptCutLevelSet
{
  public:
   xAcceptCutLevelSet(xLevelSet& ls, double _dist = 0.0, xEntityToEntity _ent = xtool::xIdentity<AOMD::mEntity*>());
   bool operator()(AOMD::mEntity* e) const;

  private:
   xLevelSet* LS;
   double dist;
   xEntityToEntity ent;
};

//! entityFilter, return true if the entity is strictly cut by the iso-zero of the level-set
struct xAcceptStrictCutLevelSet
{
  public:
   xAcceptStrictCutLevelSet(xLevelSet& ls, double _dist = 0.0, xEntityToEntity _ent = xtool::xIdentity<AOMD::mEntity*>());
   bool operator()(AOMD::mEntity* e) const;

  private:
   xLevelSet* LS;
   double dist;
   xEntityToEntity ent;
};

//! entityFilter, return true if all node of entity have levelset values less or equal to zero.
struct xAcceptStrictlyInsideLevelSet
{
  public:
   xAcceptStrictlyInsideLevelSet(xLevelSet& ls, double _dist = 0.0, xEntityToEntity _ent = xtool::xIdentity<AOMD::mEntity*>());
   bool operator()(AOMD::mEntity* e) const;

  private:
   xLevelSet* LS;
   double dist;
   xEntityToEntity ent;
};

//! entityFilter, return true if at least one  node of entity has levelset value less or equal to zero.
struct xAcceptInsideAndCrossingLevelSet
{
  public:
   xAcceptInsideAndCrossingLevelSet(xLevelSet& ls, xEntityToEntity _ent = xtool::xIdentity<AOMD::mEntity*>());
   bool operator()(AOMD::mEntity* e) const;

  private:
   xLevelSet* LS;
   xEntityToEntity ent;
};

struct xAcceptCurvatureLessThan
{
  public:
   xAcceptCurvatureLessThan(xLevelSet& ls, double _curv = 0.0, xEntityToEntity _ent = xtool::xIdentity<AOMD::mEntity*>());
   bool operator()(AOMD::mEntity* e) const;

  private:
   xLevelSet* LS;
   double curv;
   xEntityToEntity ent;
};

//! entity filter, return true if the entity is classifyed on a xDomain with name zone_name;
struct xAccept
{
  public:
   xAccept(const std::string& zone_name);
   xAccept(const int& zone_id);
   bool operator()(AOMD::mEntity* e) const;

  private:
   const int zone_id;
};

//! entity filter, return false if the entity is classifyed on a xDomain with name zone_name;
struct xAcceptNot
{
  public:
   xAcceptNot(const std::string& zone_name) : accept{zone_name} {}
   xAcceptNot(const int& zone_id) : accept{zone_id} {}
   bool operator()(AOMD::mEntity* e) const { return !accept(e); }

  private:
   xAccept accept;
};

//! entity filter, return true if the first upper adjancency of the entity is classifyed on a xDomain with name zone_name
struct xAcceptUpperAdjacency
{
  public:
   xAcceptUpperAdjacency(const std::string& zone_name) : accept{zone_name} {}
   xAcceptUpperAdjacency(const int& zone_id) : accept{zone_id} {}
   bool operator()(AOMD::mEntity* e) const;

  private:
   xAccept accept;
};

//! entity filter, return true if an entity is on a boundary of the mesh
/*! (based on, the number of entity of level n+1 attach to the current entity)
  This Can cause trouble in parallel ... Since a face can have 2 neighbor but only
  one on the current processor
*/
struct xAcceptOnBoundary
{
  public:
   bool operator()(AOMD::mEntity* e) const;
};

//! entity filter, return true if an entity is on a classified boundary condition of the mesh
/*! It return true as soon as given entity  is classified whatever its boundary condition id is.
 *  Given entity  must be of dimension less then 3 (checked in debug mode only)
 *  Top level checking would be more appropiate .?.?.
 */
struct xAcceptOnClassifiedBoundary
{
  public:
   bool operator()(AOMD::mEntity* e) const;
};

//! entity filter, return true if an entity of level 'iwhat' is on the classified boundary condition 'iclass' of the mesh
/*! All element kind are accepted.
 * If element adjacency is created for element level (AOMD adaptation or other mecahnisme) this filter return false.
 */
class xAcceptOnClassification
{
  public:
   xAcceptOnClassification(const int iclass_, int iwhat_);
   bool operator()(AOMD::mEntity* e) const;

  private:
   const int iclass, iwhat;
};

//! entity filter, return true if an entity got all its node fulfilling requirement of given predicate
/*!
 * predicate is taking  return true if x,y,z coordinate of point is in the selected volume
 */
struct xAcceptInVolume
{
  public:
   xAcceptInVolume(std::function<bool(const double& x, const double& y, const double& z)>& in_vol_);
   bool operator()(AOMD::mEntity* e) const;

  private:
   std::function<bool(const double& x, const double& y, const double& z)> in_vol;
};

//! entity filter, return true if an entity is on a boundary of a region
/*! (based on, the number of entity of level n+1 attached to the current entity, that are on the region)
 */
class xAcceptOnBoundaryOfRegion
{
  public:
   xAcceptOnBoundaryOfRegion(const xRegion& reg_);
   bool operator()(AOMD::mEntity* e) const;

  private:
   const xRegion& region;
};

//! entity filter, return true if a tag tag_ is attached on the entity and if the value associated to the tag is equal to
//! accepted_tag_val_
struct xAcceptOnIntTag
{
  public:
   xAcceptOnIntTag(const unsigned int tag_, int accepted_tag_val_);
   bool operator()(AOMD::mEntity* e) const;

  private:
   const unsigned int tag;
   const int accepted_tag_val;
};

/// A filter on mVertex * that will return true if the value of the function given at construction computed on
/*! the vertex is less than the limit given upon construction
  !*/
struct xAcceptLessThanForOneNode
{
   xAcceptLessThanForOneNode(std::function<double(AOMD::mVertex*)>& _eval, double _limit);
   bool operator()(AOMD::mEntity* e) const;

  private:
   std::function<double(AOMD::mVertex*)> eval;
   double limit;
};

// static xAcceptAll xAcceptAlltest;

}  // namespace xfem
#endif
