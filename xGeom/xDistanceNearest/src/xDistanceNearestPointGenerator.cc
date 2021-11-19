/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#include <iostream>
#include <sstream>
#include <vector>

// geom
#include "xDistanceNearestPointGenerator.h"
// xfem

#include "xLevelSet.h"
#include "xLevelSetOperators.h"
// AOMD
#include "mEntity.h"

#ifdef HAVE_CGAL  // and USE_XLEGACYSIMPLECUT
// xcut
#include "xRefCutToIsoZeroVector.h"
#endif

using namespace xfem;
using namespace AOMD;

namespace xgeom
{
#ifdef HAVE_AOMD
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xLevelSetToVectorAndTag class implementation
//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////// Constructor
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CGAL
xLevelSetToCGALAndTag::xLevelSetToCGALAndTag(xfem::xLevelSet &ls_,
                                             const xinterface::aomd::xAttachedDataManagerAOMD<char> &tag_nodes_, bool fit,
                                             double fittol)
    : refcut(xcut::xRefCutToIsoZeroVector()), ls(ls_), tagger(tag_nodes_)
{
   if (fit)
   {
      xFitToVertices fit(fittol);
      ls.accept(fit);
   }
}
xLevelSetTriToCGALEdgeAndTag::xLevelSetTriToCGALEdgeAndTag(xfem::xLevelSet &ls_,
                                                           const xinterface::aomd::xAttachedDataManagerAOMD<char> &tag_nodes_,
                                                           bool fit, double fittol)
    : xLevelSetToCGALAndTag(ls_, tag_nodes_, fit, fittol), cut_points(6)
{
}
xLevelSetTetToCGALTriAndTag::xLevelSetTetToCGALTriAndTag(xfem::xLevelSet &ls_,
                                                         const xinterface::aomd::xAttachedDataManagerAOMD<char> &tag_nodes_,
                                                         bool fit, double fittol)
    : xLevelSetToCGALAndTag(ls_, tag_nodes_, fit, fittol), cut_points(12)
{
}
xLevelSetToCGALAndTag::tagNode::tagNode(const xinterface::aomd::xAttachedDataManagerAOMD<char> &tag_nodes_)
    : tag_nodes(tag_nodes_)
{
}
#endif

/////////////////////////////////////// End constructor
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Destructor
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End Destructor
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Private methode
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CGAL
void xLevelSetToCGALAndTag::tagNode::operator()(mEntity *e) const
{
   // attachChar(e, tag_nodes, 1);
   const_cast<xinterface::aomd::xAttachedDataManagerAOMD<char> &>(tag_nodes).setData(*e) = 1;
   return;
}

void xLevelSetTriToCGALEdgeAndTag::generateEntity(EntityArg_t e, xLevelSetTriToCGALEdgeAndTag::Container_t &container) const
{
   genAndTagFromCut(*const_cast<xLevelSetTriToCGALEdgeAndTag *>(this), static_cast<mEntity *>(e), container);
   return;
}

int xLevelSetTriToCGALEdgeAndTag::cutter(std::vector<double> &vals) { return refcut.cutTriRefByLevelSet(vals, cut_points); }

void xLevelSetTriToCGALEdgeAndTag::modifyCoord(xmapping::xMapping *mapping)
{
   // due to interface virtual void eval (double u, double v, double w, double &x, double &y, double &z)
   // below implementation is possible : u,v,w arec copy of p[0],p[1],p[2] and x,y,z are reference to p[0],p[1],p[2]
   mapping->eval(cut_points[0], cut_points[1], cut_points[2], cut_points[0], cut_points[1], cut_points[2]);
   mapping->eval(cut_points[3], cut_points[4], cut_points[5], cut_points[3], cut_points[4], cut_points[5]);
   return;
}
void xLevelSetTriToCGALEdgeAndTag::modifyCoord2(xmapping::xMapping *mapping)
{
   // in 2D a tri won't give more then one iso zero edge => nothing to do
   return;
}
void xLevelSetTriToCGALEdgeAndTag::addList(xLevelSetTriToCGALEdgeAndTag::Container_t &entity_res_list)
{
   // generate CGAL geometric entity from points
   entity_res_list.emplace_back(PointRes_t(cut_points[0], cut_points[1], cut_points[2]),
                                PointRes_t(cut_points[3], cut_points[4], cut_points[5]));
   return;
}
void xLevelSetTriToCGALEdgeAndTag::addList2(xLevelSetTriToCGALEdgeAndTag::Container_t &entity_res_list)
{
   // in 2D a tri won't give more then one iso zero edge => nothing to do
   return;
}
bool xLevelSetTriToCGALEdgeAndTag::tagTouchingIso(mEntity *e, int id)
{
   if (id > 3)
   {
      // only one node is in isozero
      tagger(e->get(0, id - 4));
      return false;
   }
   else
   {
      // only one edge is in isozero
      mEntity *edge = e->get(1, id - 1);
      tagger(edge->get(0, 0));
      tagger(edge->get(0, 1));
      return true;
   }
}

void xLevelSetTetToCGALTriAndTag::generateEntity(EntityArg_t e, xLevelSetTetToCGALTriAndTag::Container_t &container) const
{
   genAndTagFromCut(*const_cast<xLevelSetTetToCGALTriAndTag *>(this), static_cast<mEntity *>(e), container);
   return;
}
int xLevelSetTetToCGALTriAndTag::cutter(std::vector<double> &vals) { return refcut.cutTetRefByLevelSet(vals, cut_points); }

void xLevelSetTetToCGALTriAndTag::modifyCoord(xmapping::xMapping *mapping)
{
   // due to interface virtual void eval (double u, double v, double w, double &x, double &y, double &z)
   // below implementation is possible : u,v,w arec copy of p[0],p[1],p[2] and x,y,z are reference to p[0],p[1],p[2]
   mapping->eval(cut_points[0], cut_points[1], cut_points[2], cut_points[0], cut_points[1], cut_points[2]);
   mapping->eval(cut_points[3], cut_points[4], cut_points[5], cut_points[3], cut_points[4], cut_points[5]);
   mapping->eval(cut_points[6], cut_points[7], cut_points[8], cut_points[6], cut_points[7], cut_points[8]);
   return;
}
void xLevelSetTetToCGALTriAndTag::modifyCoord2(xmapping::xMapping *mapping)
{
   // due to interface virtual void eval (double u, double v, double w, double &x, double &y, double &z)
   // below implementation is possible : u,v,w arec copy of p[0],p[1],p[2] and x,y,z are reference to p[0],p[1],p[2]
   mapping->eval(cut_points[9], cut_points[10], cut_points[11], cut_points[9], cut_points[10], cut_points[11]);
   return;
}
void xLevelSetTetToCGALTriAndTag::addList(xLevelSetTetToCGALTriAndTag::Container_t &entity_res_list)
{
   // generate CGAL geometric entity from points
   entity_res_list.emplace_back(PointRes_t(cut_points[0], cut_points[1], cut_points[2]),
                                PointRes_t(cut_points[3], cut_points[4], cut_points[5]),
                                PointRes_t(cut_points[6], cut_points[7], cut_points[8]));
}
void xLevelSetTetToCGALTriAndTag::addList2(xLevelSetTetToCGALTriAndTag::Container_t &entity_res_list)
{
   // generate CGAL geometric entity from points
   entity_res_list.emplace_back(PointRes_t(cut_points[0], cut_points[1], cut_points[2]),
                                PointRes_t(cut_points[6], cut_points[7], cut_points[8]),
                                PointRes_t(cut_points[9], cut_points[10], cut_points[11]));
}

bool xLevelSetTetToCGALTriAndTag::tagTouchingIso(mEntity *e, int id)
{
   if (id > 10)
   {
      // only one node is in isozero
      tagger(e->get(0, id - 11));
   }
   else if (id > 4)
   {
      // only one edge is in isozero
      mEntity *edge = e->get(1, id - 5);
      tagger(edge->get(0, 0));
      tagger(edge->get(0, 1));
   }
   else
   {
      // a face is in isozero
      mEntity *face = e->get(2, id - 1);
      tagger(face->get(0, 0));
      tagger(face->get(0, 1));
      tagger(face->get(0, 2));
      return true;
   }

   return false;
}
#endif
/////////////////////////////////////// End Private methode
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Public methode
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////// End Public methode
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

}  // namespace xgeom
