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

#ifndef _MMESH_H_
#define _MMESH_H_

#include "AOMD_SharedInfo.h"
#include "mCommonBdry.h"
#include "mEntity.h"
#include "mEntityContainer.h"
#include "mIdGenerator.h"
#include "mIterator.h"
#include "modeler.h"

#include <iosfwd>
#include <list>
#ifdef TSTT_
#include "mAttachableDataContainer.h"
#include "mEntitySet.h"
#endif

namespace AOMD
{
class mEntity;
class mFace;
class mEdge;
class mVertex;
class mMirrorVertex;
class mRegion;
class mTet;
class mHex;
class mPrism;
class AOMD_LoadBalancerCallbacks;
class AOMD_OwnerManager;
class AOMD_DataExchanger;
#ifdef TSTT_
class mEntitySet;
#endif

struct rp_int
{
   int i;
   mEntity *entity;
};

struct rp_int2
{
   int i;
   int j;
   mEntity *entity;
};

/**
 mSplitCallbacks is the base class for mesh
 refinement callbacks. mMesh class provides
 a member mMesh::refunref that allows to refine
 itself. This callbacks is used to drive mesh
 refinements and coarsenings.
*/
class mSplitCallbacks
{
  public:
   /// = 0 ok, = 1 refine, =-1 unrefine
   virtual ~mSplitCallbacks() = default;
   virtual int operator()(mEntity *) = 0;
   /// what we do when we split
   virtual void splitCallback(mEntity *) = 0;
   /// what we do when we unsplit
   virtual void unsplitCallback(mEntity *) = 0;
   void splitCallback(std::list<mEntity *> &);
   void unsplitCallback(std::list<mEntity *> &);
};

/**
 mMesh class provides an interface to
 a complete mesh that can be distributed.
*/
#ifdef TSTT_
class mMesh : public mEntitySet
#else
class mMesh
#endif
{
  protected:
   // ********************************************
   //		MESH + MODEL
   // ********************************************
   mMesh *theMeshModel;
   pGModel theSolidModel;

   // ********************************************
   //		ID
   // ********************************************
   // mesh id
   int iD;
   /// an iD generator for vertices
   mIdGenerator theIdGenerator;

   // ********************************************
   //		GENTITY RELATED
   // ********************************************
   /* flagModel	0 - null model
                   1 - GM_entityByTag
                   2 - GM_entityByID
   */

   // ********************************************
   //		CONTAINER
   // ********************************************
   /// container with 4 hash tables, all mesh entities
   mMeshEntityContainer allEntities;
   /// common boundary
   std::set<mCommonBdry *, mCommonBdryLessThanKey> allCommonBdries;

   // ********************************************
   //		ENTITY-RELATED
   // ********************************************
   mEntity *getEntity(mEntity::mType type, mVertex *v[], pGEntity g);
   unsigned int mirrorTag;

   // ********************************************
   //	NON-CONFORMING MESH ADAPTATION
   // ********************************************

   /// split and unsplit operators, should get out of this
   /// class
   //  mVertex *splitEdge(mEdge *e, mEntity *src);
   void splitTriangle(mFace *e, mEntity *src);
   mVertex *splitQuad(mFace *e, mEntity *src);
   mVertex *splitHex(mHex *e);
   void splitTet(mTet *e);
   void unsplitTriangle(mFace *e);
   void unsplitQuad(mFace *e);
   void unsplitHex(mHex *e);
   void unsplitTet(mTet *e);
   //  void unsplitEdge(mEdge *e);
   void splitgen(std::list<mEntity *> &toSplit, mSplitCallbacks &);

  public:
#ifdef TSTT_
   inline bool isMesh() { return true; }
   std::list<mEntitySet *> allEntSets;
   void addEntSet(mEntitySet *p) { allEntSets.push_back(p); }
   void removeEntSet(mEntitySet *p)
   {
      std::list<mEntitySet *>::iterator loc = std::find(allEntSets.begin(), allEntSets.end(), p);
      allEntSets.erase(loc++);
   }
#endif
   // ********************************************
   //		typedef
   // ********************************************
   /// iterators on all entities (the whole graph)
   typedef mMeshEntityContainer::iter iterall;
   /// iterators on leaves
   typedef mLeavesIterator iter;
   typedef std::set<mCommonBdry *, mCommonBdryLessThanKey>::const_iterator cbiter;

   // ********************************************
   //		CONSTRUCTOR/DESTRUCTOR
   // ********************************************
   /// constructor
   mMesh(int id = 1, pGModel model = nullptr);
   mMesh(mMesh &&) = default;
   mMesh(const mMesh &) = delete;
   mMesh &operator=(mMesh &&) = default;
   mMesh &operator=(const mMesh &) = delete;
   virtual ~mMesh();

   // ********************************************
   //		MESH LOAD RELATED
   // ********************************************
   int representationFlag;  // 0: non-one level representation
                            // 1: one-level representation
   void setRepresentationFlag(int i) { representationFlag = i; }
   int getRepresentationFlag() { return representationFlag; }

   // ********************************************
   //		MESH MODEL
   // ********************************************
   mMesh *setMeshModel(mMesh *m) { return (theMeshModel = m); }
   mMesh *getMeshModel() const { return theMeshModel; }
   pGModel getSolidModel() const { return theSolidModel; }
   void setSolidModel(pGModel m);

   // ********************************************
   //		ID
   // ********************************************
   /// we have one id generator for vertices, that function returns it
   inline const mIdGenerator &getIdGenerator() const { return theIdGenerator; }
   /// a mesh has an iD, this returns it
   inline int getId() { return iD; }
   virtual int getNewVertexId();
   virtual void setNewVertexId(int);

   // ********************************************
   //		GENTITY
   // ********************************************
   typedef pGEntity (*GEntity_FP)(pGModel, int, int);
   GEntity_FP classif_fp;
   GEntity_FP getGEntity_FP() { return classif_fp; }
   void setGEntity_FP(GEntity_FP fp) { classif_fp = fp; }
   pGEntity getGEntity(int tag, int type);
   pGEntity getGEntity(int tag, int type, pGEntity (*)(pGModel, int, int));

   // ********************************************
   //		ENTITY QUERY
   // ********************************************
   /// tells if a mesh entity exists
   virtual mEntity *find(mEntity *) const;
   /// a mesh has a container with all mesh entities
   mMeshEntityContainer &getEntities() { return allEntities; }
   /// get a vertex with id iD, returns 0 if not found
   mVertex *getVertex(int iD) const;
   /// get an edge with vertices v1 and v2, returns 0 if not found
   mEdge *getEdge(mVertex *v1, mVertex *v2) const;

   /// get an triangle with vertices v1, v2 and v3, returns 0 if not found
   mFace *getTri(mVertex *v1, mVertex *v2, mVertex *v3) const;
   /// get an quad with vertices v1, v2, v3 and v4, returns 0 if not found
   mFace *getQuad(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4) const;
   /// get a tet   with vertices v1, v2, v3 and v4, returns 0 if not found
   mTet *getTet(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4) const;
   /// get an hex with vertices v1, v2, ... and v8, returns 0 if not found
   mHex *getHex(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, mVertex *v5, mVertex *v6, mVertex *v7, mVertex *v8) const;

   //  ********************************************
   //             ONE_LEVEL related
   //  ********************************************
   /// create entities for one level representation
   mEdge *createEdge_oneLevel(mVertex *, mVertex *, pGEntity);
   mFace *createFaceWithEdges_oneLevel(mEdge *e1, mEdge *e2, mEdge *e3, pGEntity classif);
   mTet *createTetWithFaces_oneLevel(mFace *f1, mFace *f2, mFace *f3, mFace *f4, pGEntity classif);
   void DEL_updateAdj(mEntity *);
   void DEL_oneLevel(mEntity *);
   void setRepresentationToOneLevelForParallel();

   // ********************************************
   //		ENTITY CREATION
   // ********************************************
   /// add a mesh entity
   virtual void add(mEntity *);
   mEntity *copyMeshEntity(mEntity *other);

   /// create a new vertex
   mVertex *createVertex(int iD, const Trellis_Util::mPoint &p, pGEntity classif);

   /// create a new vertex
   mVertex *createVertex(int iD, const double x, const double y, const double z, pGEntity classif);

   /// create a new vertex without updating the id generator
   mVertex *createVertex_noUpdateId(int iD, const Trellis_Util::mPoint &p, pGEntity classif);

   /// create a new vertex without updating the id generator
   mVertex *createVertex_noUpdateId(int iD, const double x, const double y, const double z, pGEntity classif);

   /// create a new vertex
   mVertex *createVertex(const Trellis_Util::mPoint &p, pGEntity classif);

   /// create a new vertex
   mVertex *createVertex(const double x, const double y, const double z, pGEntity classif);

   /// create a new edge
   mEdge *createEdge(int iV1, int iV2, pGEntity classif);

   /// create a new edge
   mEdge *createEdge(mVertex *v1, mVertex *v2, pGEntity classif);

   /// create a new face
   mFace *createFaceWithVertices(int iV1, int iV2, int iV3, pGEntity classif);
   /// create a new face
   mFace *createFaceWithVertices(int iV1, int iV2, int iV3, int iV4, pGEntity classif);
   /// create a new face
   mFace *createFaceWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, pGEntity classif);
   /// create a new face
   mFace *createFaceWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, pGEntity classif);
   /// create a new face
   mFace *createFaceWithEdges(mEdge *e1, mEdge *e2, mEdge *e3, pGEntity classif, int *dir = nullptr);
   /// create a new face
   mFace *createFaceWithEdges(mEdge *e1, mEdge *e2, mEdge *e3, mEdge *e4, pGEntity classif, int *dir = nullptr);

   /// create a new tet
   mTet *createTetWithVertices(int iV1, int iV2, int iV3, int iV4, pGEntity classif);
   /// create a new tet
   mTet *createTetWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, pGEntity classif);
   /// create a new tet
   mTet *createTetWithFaces(mFace *f1, mFace *f2, mFace *f3, mFace *f4, pGEntity classif);
   /// create a new hex
   mHex *createHexWithVertices(int iV1, int iV2, int iV3, int iV4, int iV5, int iV6, int iV7, int iV8, pGEntity classif);
   /// create a new hex
   mHex *createHexWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, mVertex *v5, mVertex *v6, mVertex *v7,
                               mVertex *v8, pGEntity classif);

   /// create a new hex
   mPrism *createPrismWithVertices(int iV1, int iV2, int iV3, int iV4, int iV5, int iV6, pGEntity classif);
   /// create a new hex
   mPrism *createPrismWithVertices(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, mVertex *v5, mVertex *v6,
                                   pGEntity classif);

   // ********************************************
   //		ENTITY DELETION
   // ********************************************
   /// delete from database but does not delete the pointer
   virtual void del(mEntity *);
   /// deletes form database then deletes the pointer
   virtual void DEL(mEntity *);

   // ********************************************
   //		ENTITY STREAM
   // ********************************************
   /// reads a mesh entity on a istream
   mEntity *readMEntity(std::istream &);
   /// writes a mesh entity on a ostream
   void writeMEntity(std::ostream &, mEntity *);

   // ********************************************
   //		ADJACENCY
   // ********************************************
   void reduceToMinumumRepresentation();
   /**
      This is a very important feature of AOMD. modifyState allow
      to modify the mesh representation. We can create (delete if state = false)
      entities of dimension to starting with entities of dimension
      from using entities of dimension with

      for example modifyState (3,2) creates faces using region vertices
      for example modifyState (1,3) creates edge-region associations
      for example modifyState (0,3,false) deletes vertex-region associations
   */
   void modifyState(int from, int to, bool state = true, int with = 0);

   // ********************************************
   //		MESH CLEANUP
   // ********************************************
   void cleanup(int, mSplitCallbacks &);
   void cleanup(int);
   void cleanup_updateAdj(int);
   // this second cleanup() takes care of adjacency container
   void getEntitiesToRemove(int, std::list<mEntity *> &, std::list<mEntity *> &, std::list<mEntity *> &);

   // ********************************************
   //		MIRROR ENTITY
   // ********************************************
   /// create a mirror vertex
   // nico mMirrorVertex *createMirrorVertex ();
   /// create a mirror vertex
   // nico mMirrorVertex *createMirrorVertex (int id);
   /** for an entity e, mesh will look if there are some
       other ones that are mirrored i.e. share mirror
       vertices. l is cleared if true is returned.*/
   bool lookupForMirrorEntities(mEntity *e, std::list<mEntity *> &l);

   // ********************************************
   //		MESH INFO/ITERATOR
   // ********************************************
   /// gives the TOTAL size of entities of dimension what
   virtual int size(int what) const;

   /// iterator on all entities
   virtual iterall beginall(int what) const;
   /// end iterator on all entities
   virtual iterall endall(int what) const;
   /// iterator on leaves of the mesh (otctree leaves traversal)
   virtual iter begin(int what) const;
   /// end iterator on leaves of the mesh (otctree leaves traversal)
   virtual iter end(int what) const;
   /// iterator on entities of dimension what classified on which
   virtual mClassIterator begin(int what, int which, int onwhat) const;
   /// end iterator on entities of dimension what classified on which
   virtual mClassIterator end(int what, int which, int onwhat) const;
   /// iterator on model entities
   //    virtual giter begin() const{return allGEntities.begin();}
   /// end iterator on model entities
   //    virtual giter end  () const{return allGEntities.end();}
   /// dimension of the mesh (0,1,2 or 3).
   int getDim() const;

   /// add a mesh to another
   mMesh &operator+=(mMesh &other);

   // ********************************************
   //		COMPUTATION
   // ********************************************
   Trellis_Util::mPoint COG(mEntity *);
   Trellis_Util::mPoint closestPoint(const Trellis_Util::mPoint &pt, int dim, int tag) const;

   // ********************************************
   //	NON-CONFORMING MESH ADAPTATION
   // ********************************************
   // split/unsplit (non conformal) a list of mesh entities
   int getRefinementLevel(mEntity *e);
   int getRefinementDepth(mEntity *e);
   void refunref(mSplitCallbacks &);
   int split(std::list<mEntity *> &toSplit, mSplitCallbacks &);
   void unsplit(std::list<mEntity *> &toSplit, mSplitCallbacks &);

   // ********************************************
   //		PARALLEL
   // ********************************************
   /// a owner manager for partition boundaries
   AOMD_OwnerManager *theOwnerManager;

   // common boundary related
   mCommonBdry *getInterprocessorCommonBdry(int nbP, int *p, int dim);
   mCommonBdry *getCommonBdry(int dim, int id);
   void updateCommonBdryEntities(int what);
   /// iterator on common bdries
   virtual cbiter cbbegin() const { return allCommonBdries.begin(); }
   /// end iterator on common bdries
   virtual cbiter cbend() const { return allCommonBdries.end(); }
   int numCommonBdries()
   {
      int counter = 0;
      cbiter it = cbbegin();
      for (; it != cbend(); ++it)
      {
         counter++;
      }
      return counter;
   }
   /** Create links between partition boundaries. This is used
       in parallel as well as when we have periodic entities.
       This has to be called after all mesh modifications*/
   void bdryLinkSetup();
   // void bdryLinkSetup_oneLevel();
   void ParallelReclassifyMeshEntities(int n);
   void setCB_and_collectEntitiesToRemove(int n, std::list<mEntity *> &vtToRemove, std::list<mEntity *> &egToRemove,
                                          std::list<mEntity *> &fcToRemove);
   /** This creates a round of communications. The Data exchanger
       allows the user to easily sends and recieve datas through
       interprocessor entities or/and periodic entities. See AOMD
       example for more details
       blocking paramater : use blocking communications
       if not, non blocking (take care of dead locks then...) */
   void exchangeDataOnPartBdrys(AOMD_DataExchanger &de, bool blocking = true);

   int loadBalance(AOMD_LoadBalancerCallbacks &dm, int from, int to);
   int loadBalance_updateAdj(AOMD_LoadBalancerCallbacks &dm, int from, int to);

   /// A mesh can be also used as a model. For that, we are able to search
   /// closest points.


// ********************************************
//		ghost entity container
// ********************************************


   // ********************************************
   //		DEBUG - mDebugUtil.cc
   // ********************************************
   int numUniqueEntities(int);
   void print(int) const;
   void printAll() const;
   void printNumEntities();
   void printNumEntities_NonCfm();
// ********************************************
//		DMUM-related
// ********************************************
#ifdef DMUM
   bool DMUM_on;
   int count_trvs[4];
#endif

   // ********************************************
   //		FLEXIBLE M-REP
   // ********************************************

#if defined(FLEXDB) || defined(DMUM)
   short int MRM[4][4];
#endif

#ifdef FLEXDB
   // Mesh Representation Matrix
   // type definitions
   typedef mEdge *(*t_createE_FP)(mMesh *, mVertex *, mVertex *, pGEntity);

   typedef mFace *(*t_createTri_E_FP)(mMesh *, mEdge *, mEdge *, mEdge *, pGEntity, int *);
   typedef mFace *(*t_createQuad_E_FP)(mMesh *, mEdge *, mEdge *, mEdge *, mEdge *, pGEntity, int *);

   typedef mTet *(*t_createTet_F_FP)(mMesh *, mFace *, mFace *, mFace *, mFace *, pGEntity);
   typedef mHex *(*t_createHex_F_FP)(mMesh *, mFace *, mFace *, mFace *, mFace *, mFace *, mFace *, pGEntity);
   typedef mPrism *(*t_createPrism_F_FP)(mMesh *, mFace *, mFace *, mFace *, mFace *, mFace *, pGEntity);

   typedef mFace *(*t_createTri_V_FP)(mMesh *, mVertex *, mVertex *, mVertex *, pGEntity, int *);
   typedef mFace *(*t_createQuad_V_FP)(mMesh *, mVertex *, mVertex *, mVertex *, mVertex *, pGEntity, int *);

   typedef mTet *(*t_createTet_V_FP)(mMesh *, mVertex *, mVertex *, mVertex *, mVertex *, pGEntity);
   typedef mHex *(*t_createHex_V_FP)(mMesh *, mVertex *, mVertex *, mVertex *, mVertex *, mVertex *, mVertex *, mVertex *,
                                     mVertex *, pGEntity);
   typedef mPrism *(*t_createPrism_V_FP)(mMesh *, mVertex *, mVertex *, mVertex *, mVertex *, mVertex *, mVertex *, pGEntity);

   typedef void (*t_deleteR_FP)(mMesh *, mEntity *);
   typedef void (*t_deleteF_FP)(mMesh *, mEntity *);
   typedef void (*t_deleteE_FP)(mMesh *, mEntity *);

   // mesh mod functions
   t_createE_FP createE_FP;

   t_createTri_E_FP createTri_E_FP;
   t_createQuad_E_FP createQuad_E_FP;
   t_createTet_F_FP createTet_F_FP;
   t_createHex_F_FP createHex_F_FP;
   t_createPrism_F_FP createPrism_F_FP;

   t_createTri_V_FP createTri_V_FP;
   t_createQuad_V_FP createQuad_V_FP;
   t_createTet_V_FP createTet_V_FP;
   t_createHex_V_FP createHex_V_FP;
   t_createPrism_V_FP createPrism_V_FP;

   t_deleteR_FP deleteR_FP;
   t_deleteF_FP deleteF_FP;
   t_deleteE_FP deleteE_FP;

   void setMeshModFunctors(bool oneLevel = false);
#endif
};

#ifdef TSTT_
class mMeshLessThanKey
{
  public:
   bool operator()(mMesh *mesh1, mMesh *mesh2) const { return mesh1->getId() < mesh2->getId(); }
};
#endif
}  // namespace AOMD

#endif
