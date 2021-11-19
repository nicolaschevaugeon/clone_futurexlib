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

#ifndef _MENTITY_H_
#define _MENTITY_H_

#include <list>
#include <set>
#include <string>
#include <map>

#ifndef SIM
#include "modeler.h"
#else
#include "SimModel.h"
#endif

#include "mEntityContainer.h"
#include "mAttachableDataContainer.h"
#include "pmEntity.h"

namespace AOMD {

  class mVertex;
  class mCommonBdry;/**
     The mEntity Class is the base class for
     all mesh entities (quad, tetrahedron, hex, 
     prism, line, vertex).

     In AOMD, this entity can have any set of 
     adjacencies i.e. can have pointers to
     all other entities it is connected to or
     connects.

     An entity can be also recursively divided and
     one can also access parent and childeren
     of it. This can be used for mesh adaptation,
     multilevel methods... The mMesh class provides
     iterators on leaves of the mesh tree so that
     the user can easily access mostly refined mesh
     entities.
  */

  class mEntity : public mAttachableDataContainer
    {
#if defined(DMUM) || defined(FLEXDB)
public: 
     static mMesh* g_mesh;     
#endif


    public:
     void setPClassification(pmEntity*);
     pmEntity* getPClassification();


      // functions moved from gEntity
      double gSize();
      
      /// This is an enum of supported entities. This is useful to determine
      /// the exact type of a mesh entity. Mirror entities are used for symmetries
      /// and periodicity.
      enum mType {VERTEX,EDGE,TRI,QUAD,HEX,PRISM,PYRAMID,TET,MIRROR,FACE};
      
      enum mTopology { POINT,
                     LINE_SEGMENT,
                     POLYGON,
                     TRIANGLE, 
                     QUADRILATERAL, 
                     POLYHEDRON,
                     TETRAHEDRON, 
                     HEXAHEDRON, 
                     PRISM_, 
                     PYRAMID_, 
                     SEPTAHEDRON,
		     ALL_TOPOLOGIES,
		     UNDEFINED};
    protected:
      /// The iD is function of vertices id's, this is a
      /// hash function which is NOT UNIQUE except for vertices.
      /// Two vertices can NEVER have the same iD while 2 faces CAN.
      /// This id is NEVER consecutive and has NO interpretation
      //static unsigned int entCounter;
      int iD;
      // this is MUCH too big. This can be computed on the fly BTW
      // std::string uiD;   // string unique iD 
      /// Four adjacencies sets because our world is 3-dimensional
      /// 0 is for vertices, 1 for edges, 2 for faces and 3 for regions
      /// Note that we do not use enums for that because integers 0,1,2,3
      /// have a clear meaning here
      mAdjacencyContainer *theAdjacencies[4];
      /// Classification of the mesh entity to a model entity (model topological entity)
      pGEntity theClassification; 
      /// Classification of the mesh entity to a model entity (inter-processor entity)
      /// we use negative id's for parallel classifications (they are distingued like that from 
      /// model ones)
      mCommonBdry* theCommonBdry; 
      /// Constructor without parameters
      mEntity();
      /// This computes the mesh entity id using vertices id's
      mAttachableData *theFastData;

    public:

      /// Fast acces to data (M. Delanaye asks, and he is right !!)
      inline mAttachableData * getFastData () const { return theFastData; }
      inline void attachFastData    (mAttachableData *d) { deleteFastData(); theFastData = d; }
      inline void deleteFastData () { if(theFastData)delete theFastData; theFastData=nullptr;}
      /// virtual destructor
      virtual ~mEntity();
      /// returns the iD
      inline int getId() const {return iD;}
      inline void setId(int id) { iD = id; }
      /// returns the classification
      pGEntity getClassification() const;
      /// returns the classification
      mCommonBdry* getCommonBdry() const;
      /// Set the common boundary
      void setCommonBdry(mCommonBdry*);      /** how many sub entities of dimension "what" (for a tet *tet, 
	  tet->getNbTemplates(2) returns 4 because a tet has 4 faces). Note
	  that getNbTemplates returns always something even if the adjacency list
	  is not present.
	  getNbTemplates(what) is only valid if the dimension of the mesh entity is 
	  greater than "what" !*/
      virtual int getNbTemplates (int what) const;
      /** create the ith sub entity of level "what" using "with"
	  getTemplate always returns a new mesh entity, even if it
	  already exists in the mesh. This member is used to create
	  new sets of mesh entities.*/
      virtual mEntity* getTemplate (int ith , int what , int with) const; 
      /** gives the dimension of the entity
	  This valus is 0,1,2 or 3. Note that we do not use enums
	  here because 0,1,2 and 3 have a meaning as integers.
	  I should rename this getDim()*/
      virtual int getLevel() const = 0;
      /** gives the type of the entity (mEntity::VERTEX, mEntity::EDGE,...)
	  This of course is an enum.*/
      virtual mType getType() const = 0;
#ifdef TSTT_
      mTopology getTopo();
#endif
      /// add mEntity *m into adjacency list of the mesh entity
      void add (mEntity* m);
      /// adds e if it is not already there
      void appendUnique (mEntity* m);
      /// delete mEntity *m from adjacency list of the mesh entity
      void del (mEntity* m);
      /// lookup into adjacencies list if entity *m is present
      mEntity *find(mEntity* m)const;
      /// size(what) gives the same result as getNbTemplates for
      /// downward entities. This really checks the size of adjacency containers.
      inline int size (int what)const;
      /// iterator on adjacent entities. 
      mAdjacencyContainer::iter begin(int what);  
      mAdjacencyContainer::iter end(int what);  
      /// Random access to entities.
      /// Get the ith downward entity of level what
      inline mEntity* get(int what, int ith)const;
      /// Set classification
      void classify(pGEntity);
      /// delete all adjacency set of dimension what
      void deleteAdjacencies(int what);
      /// asks if a given adjacency list has been created
      inline int isAdjacencyCreated(int what) const;

      /// equal operator between 2 entities
      bool equal    (mEntity *) const;
      /// less than operator between 2 entities
      bool lessthan (mEntity *) const;

      void getHigherOrderUpward (int dim, mAdjacencyContainer &ents) const;
      /// compute how a mesh entity (this) use another mesh entity e
      int getUse ( mEntity *e ) const;
      /// tree functions
      /// get the leaves of a mesh entity (octree related)
      void getLeaves(std::list<mEntity*> &leaves);
      /// get all the entities under a node of the tree
      void getAllSubTree(std::list<mEntity*> &family);
      /// gat all the family : all nodes, edges, faces of the whole family
      void getAllEntitiesOfAllDimensions(std::set<mEntity*,EntityLessThanKey> &the_whole_family);
      void getAllEntitiesOfAllDimensions(std::set<mEntity*> &the_whole_family);
      void getAllEntitiesOfAllDimensions_oneLevel(std::set<mEntity*> &);
      void getFamily_oneLevel(std::list<mEntity*>&);
      void setParent(mEntity*);
      void deleteParent();
      mEntity *parent();
      mEntity *root ();
      std::string getUid () const;
      int getOwner();
      Trellis_Util::mPoint getCentroid();
      /// debug functions 
      virtual void print() const;
      /// debug functions
#ifdef HAVE_TOTALVIEW
      static int TV_ttf_display_type_ ( const AOMD::mEntity *e );
#endif
      
    };

  inline int mEntity::isAdjacencyCreated(int what) const
  {
    return (theAdjacencies[what])?1:0;
  }
  
  inline int mEntity::size (int what)const
  {    
    if(!theAdjacencies[what])
      {
	if (getLevel() > what) return getNbTemplates(what);
	else return 0;
      }
    return theAdjacencies[what]->size();
  }
    
  inline mEntity* mEntity::get(int what, int i)const
  {
    if(!theAdjacencies[what])
      return getTemplate(i,what,0);
    return (*theAdjacencies[what])[i];
  }
} // end of namespace

#ifdef HAVE_TOTALVIEW
      int TV_ttf_display_type ( const AOMD::mEntity *e );
#endif

#endif 
