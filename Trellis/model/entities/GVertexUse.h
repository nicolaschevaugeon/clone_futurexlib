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

#ifndef H_GVertexUse
#define H_GVertexUse

#include "GEntity.h"

/** A model vertex use. */
class GVertexUse : public AOMD::GEntity {
public:
  GVertexUse(SGModel *m,GVertex *v, GEdgeUsePair *eu);
  ~GVertexUse() override;

  // GEntity stuff
  GeoRep * geometry() override;
  AOMD::SBoundingBox3d bounds() const override;
  int inClosure(GEntity *ent) const override; 
  TopoType::Value type() const override;
  int dim() const override;

  /** Get the vertex this is a use of. */
  GVertex * vertex() const;

  /** Get the shell use that this vertex use is "in". */
  GShell *shell() const;

  /** Return iterator to first edge use attached to this. */
  SSListCIter<GEdgeUsePair*> firstEdgeUse() const;

  SSList<GEdge*> edges() const override;

  /** Get the edges uses that attach to this vertex use. */
  SSList<GEdgeUsePair*> edgeUses() const;

  /** Get the face uses that this vertex use is in the closure of. */
  SSList<GFaceUse*> faceUses() const;

  // these should only be called by other GEntity functions
  GVertexUse(SGModel *m, GVertex *v);
  void mergeWith(GVertexUse *);

  void replaceEdgeUse(GEdgeUsePair *old, GEdgeUsePair *newUse);
  void addEdgeUse(GEdgeUsePair *eu);
  void removeEdgeUse(GEdgeUsePair *eu);
private:
  GVertex *OwningVertex;
  SSList<GEdgeUsePair*> EdgeUses;  // wouldn't this be better to be GEdgeUse????

};

#endif
