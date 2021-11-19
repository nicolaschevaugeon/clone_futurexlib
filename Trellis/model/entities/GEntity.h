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

#ifndef H_GEntity
#define H_GEntity

#include "SModelMember.h"

#include "Range.h"
#include "SPoint3.h"

class GeoRep;
class GRegion;
class GShell;
class GFace;
class GFaceUse;
class GLoopUse;
class GEdge;
class GEdgeUsePair;
class GVertex;
class GVertexUse;
template <class T> class SSList;
class ModelEntityID;
class SString;

class SGModel;

//using AOMD::AOMD::Range;

namespace AOMD {

/** A geometric model entity. All enitites are owned by a SGModel. */
class GEntity : public SModelMember {
public:

  GEntity(SGModel *model, int tag);
  ~GEntity() override;

  /** Default name for the entity is "type tag" */
  SString name() const override;

  /** Return a renderable representation of the entity. */
  virtual GeoRep * geometry() = 0;

  /** Return the representation type of the entity. */
  virtual RepType::Value repType() const;

  /// Spatial dimension of the entity. 
  virtual int dim() const = 0;

  /** Returns true if ent is in the closure of this entity */
  virtual int inClosure(GEntity *ent) const =0; 

  /// Regions that bound this entity or that this entity bounds.
  virtual SSList<GRegion*> regions() const;

  /// Faces that bound this entity or that this entity bounds.
  virtual SSList<GFace*> faces() const;

  /// Edges that bound this entity or that this entity bounds.
  virtual SSList<GEdge*> edges() const;

  /// Vertices that bound this entity.
  virtual SSList<GVertex*> vertices() const;

  /// Underlying geometric representation of this entity.
  virtual GeomType::Value geomType() const;

  /// True if parametric space is continuous in the "dim" direction.
  virtual Logical::Value continuous(int dim) const;

  /// True if entity is periodic in the "dim" direction.
  virtual Logical::Value periodic(int dim) const;

  /// True if there are parametric degeneracies in the "dim" direction.
  virtual Logical::Value degenerate(int dim) const;

  /// Orientation of the parametric space w.r.t. the entity.
  virtual int geomDirection() const;

  /// Parametric bounds of the entity in the "i" direction.
  virtual AOMD::Range<double> parBounds(int i) const;

  /// Modeler tolerance for the entity.
  virtual double tolerance() const;

  /// True if the entity contains the given point to within tolerance.
  virtual int containsPoint(const SPoint3 &pt) const;
  /// Classify the given point w.r.t. this entity and its boundary.
  virtual int classifyPoint(const SPoint3 &pt, int chkb, GEntity **bent) const;

  /// The model owning this entity.
  SGModel *model() const;

  // should only be called by GEntity operators
  virtual GEntity * getBaseRep();
  virtual void * getNativePtr() const;
  virtual int getNativeInt() const;
};

inline SGModel * GEntity::model() const
{ return (SGModel*)SModelMember::getModel(); }



}
#endif


