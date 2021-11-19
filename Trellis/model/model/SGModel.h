/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#ifndef H_SGModel
#define H_SGModel

#include "SModel.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
#include "GRegion.h"
#include "SSList.h"
#include "SString.h"

class ModelEntityID;
class GRegion;
class GShell;
class GFace;
class GLoop;
class GEdge;
class GVertex;
class SGModel;

namespace AOMD{
class SBoundingBox3d;
}

/** A geometric model. The model is a non-manifold B-Rep. */
class SGModel : public SModel {
public:
  ~SGModel() override;

  static SGModel * createModelSMD(const SString &name);  // create from .smd file

  void stats();

  /** Returns the geometric tolerance for the entire model. */
  virtual double tolerance() const =0;

  /** Get the number of regions in this model. */
  int numRegion() const;

  /** Get the number of faces in this model. */
  int numFace() const;

  /** Get the number of edges in this model. */
  int numEdge() const;

  /** Get the number of vertices in this model. */
  int numVertex() const;

  /** Get the nth region in this model. */
  GRegion * region(int n) const;

  /** Get the nth face in this model. */
  GFace * face(int n) const;

  /** Get the nth edge in this model. */
  GEdge * edge(int n) const;

  /** Get the nth vertex in this model. */
  GVertex * vertex(int n) const;

  /** Get an iterator initialized to the first region in this model. */
  SSListCIter<GRegion*> firstRegion() const;

  /** Get an iterator initialized to the first face in this model. */
  SSListCIter<GFace*> firstFace() const;

  /** Get an iterator initialized to the first edge in this model. */
  SSListCIter<GEdge*> firstEdge() const;

  /** Get an iterator initialized to the first vertex in this model. */
  SSListCIter<GVertex*> firstVertex() const;

  /** Get the outer region for the model. At this point this always 
    returns null. */
  GRegion * outerRegion() const;

  /** Find the region with the given tag. */
  virtual GRegion * regionByTag(int n) const;
  
  /** Find the face with the given tag. */
  virtual GFace * faceByTag(int n) const;

  /** Find the edge with the given tag. */
  virtual GEdge * edgeByTag(int n) const;

  /** Find the vertex with the given tag. */
  virtual GVertex * vertexByTag(int n) const;

  virtual GRegion * regionByID(int n) const;
  virtual GFace * faceByID(int n) const;
  virtual GEdge * edgeByID(int n) const;
  virtual GVertex * vertexByID(int n) const;

  SModelMember * getMember(int type, int tag) const override;
  SModelMember * getMember(const SString &name) const override;

  void writeSMD(const SString &name);
 
  virtual void setGeomTolerance(double) {};
  void setDisplayCoherence(int ); // default is coherent

  // these should only be called by GEntity classes
  void add(GRegion *r);
  void add(GFace *f);
  void add(GEdge *e);
  void add(GVertex *v);

  void remove(GRegion *r);
  void remove(GFace *f);
  void remove(GEdge *e);
  void remove(GVertex *v);

protected:
  SGModel(const SString &name);

  virtual GVertex *createVertex(std::istream &in)=0;
  virtual GEdge *createEdge(std::istream &in)=0;
  virtual GFace *createFace(std::istream &in)=0;
  virtual GRegion *createRegion(std::istream &in)=0;

  void createOuterShell();

  // Reads in the topological information in a smd file
  void readSMD(std::istream &in);

  // data
  int d_Coherent;

  virtual void readBounds(std::istream &in);
private:

  SSList<GRegion*> RegionList;
  SSList<GFace*> FaceList;
  SSList<GEdge*> EdgeList;
  SSList<GVertex*> VertexList;

  GRegion *OuterRegion;
};


#endif


