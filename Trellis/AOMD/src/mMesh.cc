/****************************************************************************

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written
   and maintained by the Scientific Computation Research Center (SCOREC) at
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.

   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*****************************************************************************/
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>

#include "AOMD_OwnerManager.h"
#include "ParUtil.h"
#include "mAOMD.h"
#include "mBuildAdj.h"
#include "mEdge.h"
#include "mException.h"
#include "mFace.h"
#include "mHex.h"
#include "mMesh.h"
#include "mMirrorEntity.h"
#include "mPoint.h"
#include "mPrism.h"
#include "mTet.h"
#include "mVertex.h"


// added by Jie and Eunyoung, 09/2003
#include "NullEdge.h"
#include "NullFace.h"
#include "NullModel.h"
#include "NullRegion.h"
#include "NullVertex.h"

/*#include "SGModel.h"
#include "GEdge.h"
#include "GFace.h"
#include "GRegion.h"
#include "GVertex.h"
*/
#include "modeler.h"

using std::cout;
using std::endl;
using std::istream;
using std::list;
using std::ostream;
using std::set;

namespace AOMD
{
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mMesh::mMesh(int id, pGModel model)
    : iD(id), theMeshModel(nullptr), theSolidModel(model), classif_fp(GM_entityByTag), representationFlag(0)
{
   theOwnerManager = new AOMD_OwnerManager;
   mirrorTag = AOMD_Util::Instance()->lookupMeshDataId("_mirror");

#ifdef DMUM
   DMUM_on = false;
   for (int i = 0; i < 4; ++i) count_trvs[i] = 0;
#endif
}

mMesh::~mMesh()
{
   mEntity *tmp;
   for (int i = 3; i >= 0; i--)
   {
      iterall it = beginall(i);
      iterall ite = endall(i);
      for (; it != ite;)
      {
         tmp = *it;
         ++it;
         delete tmp;
      }
   }
   delete theOwnerManager;
#ifdef TSTT_
   allEntSets.clear();
#endif
}
int mMesh::size(int what) const { return allEntities.size(what); }

mMesh::iterall mMesh::beginall(int what) const { return allEntities.begin(what); }

mMesh::iterall mMesh::endall(int what) const { return allEntities.end(what); }

mMesh::iter mMesh::begin(int what) const { return mLeavesIterator(beginall(what), endall(what)); }

mMesh::iter mMesh::end(int what) const { return mLeavesIterator(endall(what), endall(what)); }

mClassIterator mMesh::begin(int what, int which, int onwhat) const
{
   return mClassIterator(beginall(what), endall(what), which, onwhat);
}

mClassIterator mMesh::end(int what, int which, int onwhat) const
{
   return mClassIterator(endall(what), endall(what), which, onwhat);
}

mVertex *mMesh::getVertex(int theId) const
{
   static mVertex v(0, Trellis_Util::mPoint(0, 0, 0), nullptr);
   v.setIdNoRand(theId);
   return (mVertex *)find(&v);
}

mEdge *mMesh::getEdge(mVertex *v1, mVertex *v2) const
{
   static mEdge e(nullptr, nullptr, nullptr);
   e.setVertices(v1, v2);
   return (mEdge *)find(&e);
}

mFace *mMesh::getTri(mVertex *v1, mVertex *v2, mVertex *v3) const
{
   mFace f(v1, v2, v3, nullptr);
   return (mFace *)find(&f);
}

mFace *mMesh::getQuad(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4) const
{
   mFace f(v1, v2, v3, v4, nullptr);
   return (mFace *)find(&f);
}

mTet *mMesh::getTet(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4) const
{
   mTet t(v1, v2, v3, v4, nullptr);
   return (mTet *)find(&t);
}

mHex *mMesh::getHex(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, mVertex *v5, mVertex *v6, mVertex *v7, mVertex *v8) const
{
   mHex h(v1, v2, v3, v4, v5, v6, v7, v8, nullptr);
   return (mHex *)find(&h);
}

mVertex *mMesh::createVertex(int theId, const double x, const double y, const double z, pGEntity classif)
{
   return createVertex(theId, Trellis_Util::mPoint(x, y, z), classif);
}

mVertex *mMesh::createVertex(int theId, const Trellis_Util::mPoint &p, pGEntity classif)
{
   mVertex *theNewVertex = getVertex(theId);

   if (theNewVertex)
   {
      //	ParUtil::Instance()->Msg(ParUtil::WARNING,"creating an existing vertex %d\n",theId);
      std::cout << "P" << ParUtil::Instance()->rank() << ": WARNING: creating an existing vertex " << theId << std::endl;
      return theNewVertex;
   }

   theNewVertex = new mVertex(theId, p, classif);
   //    theNewVertex->print();
   theIdGenerator.setMaxValue(theId);

   setNewVertexId(theId);
   add(theNewVertex);
   return theNewVertex;
}
mVertex *mMesh::createVertex_noUpdateId(int theId, const Trellis_Util::mPoint &p, pGEntity classif)
{
   mVertex *theNewVertex = new mVertex(theId, p, classif);
   add(theNewVertex);
   return theNewVertex;
}

mVertex *mMesh::createVertex_noUpdateId(int theId, const double x, const double y, const double z, pGEntity classif)
{
   return createVertex_noUpdateId(theId, Trellis_Util::mPoint(x, y, z), classif);
}

bool mMesh::lookupForMirrorEntities(mEntity *e, std::list<mEntity *> &mirrors)
{
   mVertex *v;
   mMirrorVertex *mv[256];
   mAttachableMirrorVertex *att;

   if (e->getLevel() == 0)
   {
      v = (mVertex *)e;
      if (!(att = (mAttachableMirrorVertex *)v->getData(mirrorTag)))
         return false;
      else
      {
         mv[0] = (mMirrorVertex *)att->e;
      }
      for (int i = 0; i < mv[0]->nbCopies(); i++)
      {
         mVertex *v1 = mv[0]->getCopy(i);
         if (v1 != v) mirrors.push_back(v1);
      }
      return true;
   }

   for (int i = 0; i < e->size(0); i++)
   {
      v = (mVertex *)e->get(0, i);
      //      if(!(att = (mAttachableEntity*)v->getData("mirror")))return false;
      if (!(att = (mAttachableMirrorVertex *)v->getData(mirrorTag)))
         return false;
      else
      {
         mv[i] = (mMirrorVertex *)att->e;
      }
   }

   mirrors.clear();
   if (e->getType() == mEntity::EDGE)
   {
      for (int i = 0; i < mv[0]->nbCopies(); i++)
      {
         mVertex *v1 = mv[0]->getCopy(i);
         for (int j = 0; j < mv[1]->nbCopies(); j++)
         {
            mVertex *v2 = mv[1]->getCopy(j);
            mEdge *med = getEdge(v1, v2);
            if (med && med != e) mirrors.push_back(med);
         }
      }
   }
   else
   {
      throw 1;
   }
   return true;
}

int mMesh ::getNewVertexId() { return theIdGenerator.generateId(); }

void mMesh ::setNewVertexId(int theId) { theIdGenerator.setMaxValue(theId); }



mVertex *mMesh::createVertex(const double x, const double y, const double z, pGEntity classif)
{
   return createVertex(getNewVertexId(), Trellis_Util::mPoint(x, y, z), classif);
}
mVertex *mMesh::createVertex(const Trellis_Util::mPoint &p, pGEntity classif)
{
   return createVertex(getNewVertexId(), p, classif);
}

mEdge *mMesh::createEdge(int i1, int i2, pGEntity classif)
{
   mVertex *v1 = getVertex(i1);
   mVertex *v2 = getVertex(i2);
   if (!v1 || !v2)
   {
      char text[256];
      sprintf(text, "unknown verex id's %d(%p) %d(%p) in new edge creation\n", i1, v1, i2, v2);
      throw new mException(__LINE__, __FILE__, text);
   }
   return createEdge(v1, v2, classif);
}

void mMesh ::DEL(mEntity *e)
{
   del(e);
   delete e;
}

mEdge *mMesh::createEdge(mVertex *v1, mVertex *v2, pGEntity classif)
{
   mEdge *theNewEdge = new mEdge(v1, v2, classif);
   //        theNewEdge->print();
   add(theNewEdge);
   return theNewEdge;
}

mFace *mMesh::createFaceWithVertices(int i1, int i2, int i3, pGEntity classif)
{
   mVertex *v1 = getVertex(i1);
   mVertex *v2 = getVertex(i2);
   mVertex *v3 = getVertex(i3);
   if (!v1 || !v2 || !v3)
   {
      throw new mException(__LINE__, __FILE__, "unknown verex id's in new face creation");
   }
   return createFaceWithVertices(v1, v2, v3, classif);
}

mFace *mMesh::createFaceWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, pGEntity classif)
{
   mFace *theNewFace = new mFace(v1, v2, v3, classif);
   //    theNewFace->print();
   add(theNewFace);
   return theNewFace;
}

mFace *mMesh::createFaceWithVertices(int i1, int i2, int i3, int i4, pGEntity classif)
{
   mVertex *v1 = getVertex(i1);
   mVertex *v2 = getVertex(i2);
   mVertex *v3 = getVertex(i3);
   mVertex *v4 = getVertex(i4);
   if (!v1 || !v2 || !v3 || !v4)
   {
      char text[256];
      sprintf(text, "trying to create a quad with vertices %d %d %d %d, search gives %p %p %p %p\n", i1, i2, i3, i4, v1, v2, v3,
              v4);
      throw new mException(__LINE__, __FILE__, text);
   }

   return createFaceWithVertices(v1, v2, v3, v4, classif);
}

mFace *mMesh::createFaceWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, pGEntity classif)
{
   if (!v1 || !v2 || !v3 || !v4)
   {
      char text[256];
      sprintf(text, "trying to create a quad with %p %p %p %p\n", v1, v2, v3, v4);
      throw new mException(__LINE__, __FILE__, text);
   }

   mFace *theNewFace = new mFace(v1, v2, v3, v4, classif);
   //    theNewFace->print();
   add(theNewFace);
   return theNewFace;
}

mTet *mMesh::createTetWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, pGEntity classif)
{
   mTet *theNewTet;
   theNewTet = new mTet(v1, v2, v3, v4, classif);
   //    theNewTet->print();
   add(theNewTet);
   return theNewTet;
}

mTet *mMesh::createTetWithVertices(int i1, int i2, int i3, int i4, pGEntity classif)
{
   mVertex *v1 = getVertex(i1);
   mVertex *v2 = getVertex(i2);
   mVertex *v3 = getVertex(i3);
   mVertex *v4 = getVertex(i4);
   if (!v1 || !v2 || !v3 || !v4)
   {
      throw new mException(__LINE__, __FILE__, "unknown verex id's in new face creation");
   }
   return createTetWithVertices(v1, v2, v3, v4, classif);
}

mTet *mMesh::createTetWithFaces(mFace *f1, mFace *f2, mFace *f3, mFace *f4, pGEntity classif)
{
   mTet *theNewTet;
   theNewTet = new mTet(f1, f2, f3, f4, classif);
   //    theNewTet->print();
   add(theNewTet);
   return theNewTet;
}

mFace *mMesh::createFaceWithEdges(mEdge *e1, mEdge *e2, mEdge *e3, mEdge *e4, pGEntity classif, int *dir)
{
   mFace *theNewFace;
   theNewFace = new mFace(e1, e2, e3, e4, classif, dir);
   //    theNewFace->print();
   add(theNewFace);
   return theNewFace;
}

mFace *mMesh::createFaceWithEdges(mEdge *e1, mEdge *e2, mEdge *e3, pGEntity classif, int *dir)
{
   mFace *theNewFace;
   theNewFace = new mFace(e1, e2, e3, classif, dir);
   //    theNewFace->print();
   add(theNewFace);
   return theNewFace;
}

mHex *mMesh::createHexWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, mVertex *v5, mVertex *v6, mVertex *v7,
                                   mVertex *v8, pGEntity classif)
{
   mHex *theNewHex = new mHex(v1, v2, v3, v4, v5, v6, v7, v8, classif);
   //    theNewHex->print();
   add(theNewHex);
   return theNewHex;
}

mPrism *mMesh::createPrismWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, mVertex *v5, mVertex *v6,
                                       pGEntity classif)
{
   mPrism *theNewPrism = new mPrism(v1, v2, v3, v4, v5, v6, classif);
   //    theNewPrism->print();
   add(theNewPrism);
   return theNewPrism;
}

mHex *mMesh::createHexWithVertices(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, pGEntity classif)
{
   mVertex *v1 = getVertex(i1);
   mVertex *v2 = getVertex(i2);
   mVertex *v3 = getVertex(i3);
   mVertex *v4 = getVertex(i4);
   mVertex *v5 = getVertex(i5);
   mVertex *v6 = getVertex(i6);
   mVertex *v7 = getVertex(i7);
   mVertex *v8 = getVertex(i8);
   if (!v1 || !v2 || !v3 || !v4 || !v5 || !v6 || !v7 || !v8)
   {
      throw new mException(__LINE__, __FILE__, "unknown verex id's in new hex creation");
   }
#ifdef _DEBUG_
#endif
   return createHexWithVertices(v1, v2, v3, v4, v5, v6, v7, v8, classif);
}

mPrism *mMesh::createPrismWithVertices(int i1, int i2, int i3, int i4, int i5, int i6, pGEntity classif)
{
   mVertex *v1 = getVertex(i1);
   mVertex *v2 = getVertex(i2);
   mVertex *v3 = getVertex(i3);
   mVertex *v4 = getVertex(i4);
   mVertex *v5 = getVertex(i5);
   mVertex *v6 = getVertex(i6);
   if (!v1 || !v2 || !v3 || !v4 || !v5 || !v6)
   {
      throw new mException(__LINE__, __FILE__, "unknown verex id's in new prism creation");
   }
   return createPrismWithVertices(v1, v2, v3, v4, v5, v6, classif);
}

void mMesh::modifyState(int from, int to, bool state, int with)
{
   if (!state)
   {
      std::for_each(beginall(from), endall(from), deleteAdjFunctor(to));
   }
   else if (from > to)
   {
      std::for_each(beginall(from), endall(from), createDownwardFunctor(to, with, this));
   }
   else if (from < to)
   {
      std::for_each(begin(to), end(to), createUpwardFunctor(from));
   }
}

pGEntity mMesh::getGEntity(int tag, int type) { return getGEntity(tag, type, classif_fp); }

pGEntity mMesh::getGEntity(int tag, int type, pGEntity (*pf)(pGModel, int, int))
{
   pGEntity theSolidModelEntity;
   bool flag_NullModel = false;
   if (!theSolidModel) theSolidModel = new NullModel();
   if (!strcmp((char *)(theSolidModel->modeler()), "null")) flag_NullModel = true;
   if (!flag_NullModel)
   {
      theSolidModelEntity = pf(theSolidModel, type, tag);
      if (!theSolidModelEntity)
      {
         printf("model entity with dim %d and tag %d not found...\n", type, tag);
         throw;
      }
   }
   else
   {
      theSolidModelEntity = GM_entityByTag(theSolidModel, type, tag);
      if (!theSolidModelEntity)
      {
         // if the model entity has not been defined, define it
         switch (type)
         {
            case 0:  // model vertex
               theSolidModelEntity = new NullVertex((NullModel *)theSolidModel, tag);
               break;
            case 1:  // model edge
               theSolidModelEntity = new NullEdge((NullModel *)theSolidModel, tag);
               break;
            case 2:  // model face
               theSolidModelEntity = new NullFace((NullModel *)theSolidModel, tag);
               break;
            case 3:  // model region
               theSolidModelEntity = new NullRegion((NullModel *)theSolidModel, tag);
               break;
         }
      }
   }
   return theSolidModelEntity;
}

mCommonBdry *mMesh::getCommonBdry(int id, int l)
{
   //    int tag=-1;
   mCommonBdry gg(id, l);
   std::set<mCommonBdry *, mCommonBdryLessThanKey>::const_iterator it = allCommonBdries.find(&gg);
   if (it != allCommonBdries.end()) return *it;
   mCommonBdry *pg = new mCommonBdry(id, l);
   allCommonBdries.insert(pg);

   return pg;
}

mEntity *mMesh::find(mEntity *e) const { return allEntities.find(e); }

void mMesh::add(mEntity *e) { allEntities.add(e); }

void mMesh::del(mEntity *e) { allEntities.del(e); }

mEntity *mMesh::getEntity(mEntity::mType type, mVertex *v[], pGEntity g)
{
   mEntity *e;
   switch (type)
   {
      case mEntity::HEX:
         e = (mEntity *)getHex(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
         if (!e) e = (mEntity *)createHexWithVertices(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], g);
         break;
      case mEntity::TET:
         e = (mEntity *)getTet(v[0], v[1], v[2], v[3]);
         if (!e) e = (mEntity *)createTetWithVertices(v[0], v[1], v[2], v[3], g);
         break;
      case mEntity::QUAD:
         e = (mEntity *)getQuad(v[0], v[1], v[2], v[3]);
         if (!e) e = (mEntity *)createFaceWithVertices(v[0], v[1], v[2], v[3], g);
         break;
      case mEntity::TRI:
         e = (mEntity *)getTri(v[0], v[1], v[2]);
         if (!e) e = (mEntity *)createFaceWithVertices(v[0], v[1], v[2], g);
         break;
      case mEntity::EDGE:
         e = (mEntity *)getEdge(v[0], v[1]);
         if (!e) e = (mEntity *)createEdge(v[0], v[1], g);
         break;
   }
   return e;
}

mEntity *mMesh::copyMeshEntity(mEntity *other)
{
   mVertex *v[12];
   pGEntity g;

   if (!other->getClassification())
   {
      ParUtil::Instance()->Msg(ParUtil::WARNING, "Unclassified entity %d %d\n", other->getLevel(), other->getId());
      g = getGEntity(1000000, other->getLevel());
   }
   else
   {
      g = getGEntity(GEN_tag(other->getClassification()), GEN_type(other->getClassification()));
   }

   if (other->getLevel() == 0)
   {
      if (mVertex *v = getVertex(other->getId()))
      {
         return (mEntity *)v;
      }
      else
      {
         mVertex *vv = (mVertex *)other;
         Trellis_Util::mPoint p(vv->point());
         return (mEntity *)createVertex(vv->getId(), p, g);
      }
   }

   for (int i = 0; i < other->size(0); i++)
   {
      v[i] = getVertex(other->get(0, i)->getId());
      if (!v[i]) ParUtil::Instance()->Msg(ParUtil::ERROR, "Unable to find vertex %d\n", other->get(0, i)->getId());
   }
   mEntity *e = nullptr;
   e = getEntity(other->getType(), v, g);
   return e;
}

mMesh &mMesh::operator+=(mMesh &other)
{
   for (int i = 0; i < 4; i++)
   {
      for (iterall it = other.beginall(i); it != other.endall(i); ++it)
      {
         mEntity *e = copyMeshEntity(*it);
         if (!find(e)) add(e);
      }
   }
   return *this;
}

void mMesh::writeMEntity(ostream &o, mEntity *m)
{
   o << m->getType() << " ";
   for (int i = 0; i < m->getNbTemplates(0); i++)
   {
      o << m->get(0, i)->getId() << " ";
   }
   o << "\n";
}

mEntity *mMesh::readMEntity(istream &is)
{
   int typ, iv[10];
   is >> typ;
   // printf("reading a mesh entity of typ %d\n",typ);
   switch (typ)
   {
      case mEntity::EDGE:
         is >> iv[0] >> iv[1];
         //      if(ParUtil::Instance()->master())
         //	printf("it's an edge (%d %d)\n",iv[0],iv[1]);
         return getEdge(getVertex(iv[0]), getVertex(iv[1]));
      case mEntity::TRI:
         is >> iv[0] >> iv[1] >> iv[2];
         return getTri(getVertex(iv[0]), getVertex(iv[1]), getVertex(iv[2]));
      case mEntity::QUAD:
         is >> iv[0] >> iv[1] >> iv[2] >> iv[3];
         return getQuad(getVertex(iv[0]), getVertex(iv[1]), getVertex(iv[2]), getVertex(iv[3]));
      case mEntity::TET:
         is >> iv[0] >> iv[1] >> iv[2] >> iv[3];
         return getTet(getVertex(iv[0]), getVertex(iv[1]), getVertex(iv[2]), getVertex(iv[3]));
      case mEntity::HEX:
         is >> iv[0] >> iv[1] >> iv[2] >> iv[3] >> iv[4] >> iv[5] >> iv[6] >> iv[7];
         return getHex(getVertex(iv[0]), getVertex(iv[1]), getVertex(iv[2]), getVertex(iv[3]), getVertex(iv[4]), getVertex(iv[5]),
                       getVertex(iv[6]), getVertex(iv[7]));
      default:
       std::cout << "Warning mType "<<  typ << " not treated in switch "<< __FILE__ << " " << __LINE__ << std::endl;
       return nullptr;
   }
}

/*
insure that all inter-processor entities have
the right classification
*/

void mMesh::updateCommonBdryEntities(int what)
{
   for (iterall it = beginall(what); it != endall(what); ++it)
   {
      mEntity *e = *it;
      mCommonBdry *g = e->getCommonBdry();
      if (g)
      {
         for (int dim = 0; dim < what; dim++)
         {
            if (e->isAdjacencyCreated(dim))
            {
               for (int j = 0; j < e->size(dim); j++)
               {
                  mEntity *sub = e->get(dim, j);
                  if (!sub->getCommonBdry()) sub->setCommonBdry(g);
               }
            }
         }
      }
   }
}

void mMesh::reduceToMinumumRepresentation()
{
   int dim = (size(3)) ? 3 : 2;
   for (int DIM = dim - 1; DIM > 0; DIM--)
   {
      std::list<mEntity *> todel;
      iterall it = beginall(DIM);
      iterall itend = endall(DIM);
      for (; it != itend; ++it)
      {
         if (GEN_type((*it)->getClassification()) != (DIM) && !(*it)->getCommonBdry()) todel.push_back(*it);
      }
      for (std::list<mEntity *>::const_iterator itl = todel.begin(); itl != todel.end(); ++itl) DEL(*itl);
   }
}

int mMesh::getDim() const
{
   if (size(3)) return 3;
   if (size(2)) return 2;
   if (size(1)) return 1;
   return 0;
}

Trellis_Util::mPoint AOMD::mMesh::closestPoint(const Trellis_Util::mPoint &pt, int dim, int tag) const
{
   mClassIterator it = begin(dim, tag, dim);
   mClassIterator ite = end(dim, tag, dim);

   if (it == ite) throw;

   double mind = 1.e20;

   Trellis_Util::mPoint pp;

   /// model edge
   if (dim == 1)
   {
      while (it != ite)
      {
         mEdge *e = (mEdge *)*it;
         Trellis_Util::mVector v1(e->vertex(0)->point(), e->vertex(1)->point());
         Trellis_Util::mVector v2(e->vertex(0)->point(), pt);
         double nv1 = v1.normValue();
         double nv2 = v2.normValue();
         Trellis_Util::mVector v1Xv2(v1 % v2);
         double d = fabs(v1Xv2.normValue() * nv2);
         if (d < mind)
         {
            double xx = fabs(v1 * v2);
            mind = d;
            pp = e->vertex(0)->point() * (1. - xx) + e->vertex(1)->point() * xx;
         }
         ++it;
      }
   }
   else
   {
      throw;
   }
   return pp;
}

void mMesh::setSolidModel(pGModel m) { theSolidModel = m; }

}  // namespace AOMD
