#ifndef XAOMDENTITYUTIL_H
#define XAOMDENTITYUTIL_H
#include <cassert>
#include <iostream>
// aomd trellis
#include "AOMD_Internals.h"
#include "mEntity.h"
#include "mVertex.h"
// xtensor
#include "xPoint.h"
// xmapping
#include "xReferenceElement.h"

namespace xinterface
{
namespace aomd
{
/// return the constitutive entity of entity e
/*! d is the level of constitutive entity needed
  getSubentities add things to l
  So make sure it is empty before calling it
  It's done on purpose so don't change this behaviour without changing at least xcut::xPhysSurfByTagging::construct_ */
inline void getSubentities(AOMD::mEntity *e, const int d, std::vector<AOMD::mEntity *> &l)
{
   // if the constitutive entity searched are of same level as e, or of higher level
   // this call don't have much sens
   // return imediatly leaving the caller with it's vector
   // => caller have the responsability of checking the corectness of l
   if (d >= e->getLevel()) return;
   // number of Constitutive entity of entity e
   const int n = e->size(d);
   // loop on constitutive entity
   for (int i = 0; i < n; ++i) l.push_back(e->get(d, i));
   return;
}

/// return the list of entity on the boundary of the support the entity in.
/*! Only works for 2d mesh of simplexes. */
inline void getSupportBoundary(AOMD::mEntity *in, std::set<AOMD::mEntity *> &supportBoundary)
{
   supportBoundary.clear();
   int nface = in->size(2);
   if (in->getLevel() == 0)
   {
      for (int i = 0; i < nface; ++i)
      {
         AOMD::mEntity *currentface = in->get(2, i);
         for (int k = 0; k < 3; ++k)
         {
            AOMD::mEntity *currentedge = currentface->get(1, k);
            if (currentedge->size(2) == 1)
               supportBoundary.insert(currentedge);
            else
            {
               AOMD::mEntity *v0 = currentedge->get(0, 0);
               AOMD::mEntity *v1 = currentedge->get(0, 1);
               if (v0 != in && v1 != in) supportBoundary.insert(currentedge);
            }
         }
      }
   }
   else if (in->getLevel() == 1)
   {
      if (nface == 1) supportBoundary.insert(in);
      for (int i = 0; i < nface; ++i)
      {
         AOMD::mEntity *currentface = in->get(2, i);
         for (int k = 0; k < 3; ++k)
         {
            AOMD::mEntity *currentedge = currentface->get(1, k);
            if (currentedge != in) supportBoundary.insert(currentedge);
         }
      }
   }
   else
   {
      std::cout << "Warning : not  a node  nor an edge in get support boundary " << std::endl;
   }
}

inline xtensor::xPoint centroid(const AOMD::mEntity &e)
{
   if (e.getLevel() == 0) return static_cast<const AOMD::mVertex &>(e).point();
   xtensor::xPoint g;
   for (int i = 0; i < e.size(0); ++i) g += static_cast<const AOMD::mVertex *>(e.get(0, i))->point();
   return g *= 1. / e.size(0);
}

//! Return an std::vector containing the vertices of entity e.
inline std::vector<const AOMD::mVertex *> getVertices(const AOMD::mEntity &e)
{
   if (e.getLevel() == 0) return std::vector<const AOMD::mVertex *>{static_cast<const AOMD::mVertex *>(&e)};
   const size_t nv = size_t(e.size(0));
   std::vector<const AOMD::mVertex *> vertices(nv);
   for (size_t i = 0; i < nv; ++i) vertices[i] = static_cast<const AOMD::mVertex *>(e.get(0, int(i)));
   return vertices;
}

//! Return an std::vector containing the point associated to each vertex of entity e.
inline std::vector<xtensor::xPoint> getPoints(const AOMD::mEntity &e)
{
   if (e.getLevel() == 0) return std::vector<xtensor::xPoint>{static_cast<const AOMD::mVertex &>(e).point()};
   const size_t nv = size_t(e.size(0));
   std::vector<xtensor::xPoint> points(nv);
   for (size_t i = 0; i < nv; ++i) points[i] = static_cast<AOMD::mVertex *>(e.get(0, int(i)))->point();
   return points;
}
//! Return an  xPoint associated to vertex e.
inline xtensor::xPoint getPoint(const AOMD::mEntity &e)
{
   assert(e.getLevel() == 0);
   return static_cast<xtensor::xPoint>(static_cast<const AOMD::mVertex &>(e).point());
}

inline mType mapTomType(xmapping::xReferenceElementType in)
{
   using namespace xmapping;
   switch (in)
   {
      case xmapping::xReferenceElementType::VERTEX:
         return VERTEX;
      case xmapping::xReferenceElementType::EDGE:
         return EDGE;
      case xmapping::xReferenceElementType::TRI:
         return TRI;
      case xmapping::xReferenceElementType::QUAD:
         return QUAD;
      case xmapping::xReferenceElementType::HEX:
         return HEX;
      case xmapping::xReferenceElementType::PRISM:
         return PRISM;
      case xmapping::xReferenceElementType::PYRAMID:
         return PYRAMID;
      case xmapping::xReferenceElementType::TET:
         return TET;
         // Might be usefull if we insert new xReferenceElementType that are not mappable to Trellis mType :
      default:
         throw -5456;
   }
   throw;
}

inline xmapping::xReferenceElementType mapToxReferenceElementType(AOMD::mEntity::mType in)
{
   using namespace xmapping;
   switch (in)
   {
      case AOMD::mEntity::mType::VERTEX:
         return xmapping::xReferenceElementType::VERTEX;
      case AOMD::mEntity::mType::EDGE:
         return xmapping::xReferenceElementType::EDGE;
      case AOMD::mEntity::mType::TRI:
         return xmapping::xReferenceElementType::TRI;
      case AOMD::mEntity::mType::QUAD:
         return xmapping::xReferenceElementType::QUAD;
      case AOMD::mEntity::mType::HEX:
         return xmapping::xReferenceElementType::HEX;
      case AOMD::mEntity::mType::PRISM:
         return xmapping::xReferenceElementType::PRISM;
      case AOMD::mEntity::mType::PYRAMID:
         return xmapping::xReferenceElementType::PYRAMID;
      case AOMD::mEntity::mType::TET:
         return xmapping::xReferenceElementType::TET;
         // Might be usefull if we insert new mType are introduced that are not mapable to xReferenceElementType :
      default:
         throw -684547;
   }
   throw;
}

inline xmapping::xReferenceElementType getRefElementType(const AOMD::mEntity &e)
{
   // AOMD::mEntity::mType mtype = M_GetElementType(const_cast<AOMD::mEntity *>(&e));
   AOMD::mEntity::mType mtype = (const_cast<AOMD::mEntity *>(&e))->getType();
   return mapToxReferenceElementType(mtype);
}

/*! A crud function which print in given stream all topological information of an AOMD mEntity using node ID.
 * If deep argument is true all declared adjacency of the given entity is also printed by using a recursive call
 * to this function.
 * Adress and ID of given entity are also printed
 */
inline void printEntity(std::ostream &os, const AOMD::mEntity *e, bool deep)
{
   if (deep) os << std::endl;
   int n;
   switch (e->getLevel())
   {
      case 0:
      {
         os << "node " << e << " " << e->getId() << std::endl;

         for (int i = 1; deep && i < 4; ++i)
         {
            if (e->isAdjacencyCreated(i))
            {
               os << " adjancy " << i << std::endl;
               for (int k = 0, m = e->size(i); k < m; ++k) printEntity(os, e->get(i, k), false);
            }
         }
         break;
      }
      case 1:
      {
         n = 2;
         os << "edge " << e << " " << e->getId();
         if (e->isAdjacencyCreated(0))
         {
            if (e->size(0) == n)
            {
               os << " (" << e->get(0, 0)->getId();
               for (int k = 1; k < n; ++k) os << "-" << e->get(0, k)->getId();
               os << ")" << std::endl;
            }
            else
               os << " (wrong number of nodes !! sould be " << n << " but is " << e->size(0) << ")" << std::endl;
         }

         for (int i = 1; deep && i < 4; ++i)
         {
            if (e->isAdjacencyCreated(i))
            {
               os << " adjancy " << i << std::endl;
               for (int k = 0, m = e->size(i); k < m; ++k) printEntity(os, e->get(i, k), false);
            }
         }
         break;
      }
      case 2:
      {
         n = 3;
         os << "tria " << e << " " << e->getId();
         if (e->isAdjacencyCreated(0))
         {
            if (e->size(0) == n)
            {
               os << " (" << e->get(0, 0)->getId();
               for (int k = 1; k < n; ++k) os << "-" << e->get(0, k)->getId();
               os << ")" << std::endl;
            }
            else
               os << " (wrong number of nodes !! sould be " << n << " but is " << e->size(0) << ")" << std::endl;
         }

         for (int i = 1; deep && i < 4; ++i)
         {
            if (e->isAdjacencyCreated(i))
            {
               os << " adjancy " << i << std::endl;
               for (int k = 0, m = e->size(i); k < m; ++k) printEntity(os, e->get(i, k), false);
            }
         }
         break;
      }
      case 3:
      {
         n = 4;
         os << "tet " << e << " " << e->getId();
         if (e->isAdjacencyCreated(0))
         {
            if (e->size(0) == n)
            {
               os << " (" << e->get(0, 0)->getId();
               for (int k = 1; k < n; ++k) os << "-" << e->get(0, k)->getId();
               os << ")" << std::endl;
            }
            else
               os << " (wrong number of nodes !! sould be " << n << " but is " << e->size(0) << ")" << std::endl;
         }

         for (int i = 1; deep && i < 4; ++i)
         {
            if (e->isAdjacencyCreated(i))
            {
               os << " adjancy " << i << std::endl;
               for (int k = 0, m = e->size(i); k < m; ++k) printEntity(os, e->get(i, k), false);
            }
         }
         break;
      }
   }
}
}  // end namespace aomd
}  // namespace xinterface

#endif  // XAOMDENTITYUTIL_H
