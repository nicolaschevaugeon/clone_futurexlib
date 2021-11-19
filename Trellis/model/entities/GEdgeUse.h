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
#ifndef H_GEdgeUse
#define H_GEdgeUse

#include "GEntity.h"

class GLoopUse;
class GEdge;
class GVertex;
class GVertexUse;
class GEdgeUsePair;

/** One side of an edge use. In the full radial edge data structure this
  would be an edge use. 
  */
class GEdgeUse : public AOMD::GEntity {
public:
  GEdgeUse(SGModel *m, GEdgeUsePair *ep, GLoopUse *l, GVertex *v0);
  ~GEdgeUse() override;

  // GEntity stuff
  GeoRep * geometry() override;
  AOMD::SBoundingBox3d bounds() const override;
  int inClosure(GEntity *ent) const override; // if ent is in closure of this 
  TopoType::Value type() const override;
  int dim() const override;

  /** Get the edge. */
  GEdge *edge() const;
  
  /** Get the shell use. */
  GShell *shell() const;

  /** Get the loop use. */
  GLoopUse *loopUse() const;

  /** Get the vertex use. */
  GVertexUse *vertexUse() const;

  /** Get the direction of this edge use side w.r.t the edge. */ 
  int dir() const;

  /** Get the owning edge use pair. */
  GEdgeUsePair *use() const;

  /** Get the mate to this side. */
  GEdgeUse *otherSide();

  //
  void setLoopUse(GLoopUse *lu);
  void setVertexUse(GVertexUse *vu);

  GEdgeUse(SGModel *m, GEdgeUsePair *e);

  // debugging
  void info() const;
protected:

private:
  GEdgeUsePair * Use;
  GLoopUse * LoopUse; 
  GVertexUse *VertexUse;  // vertex use at the head of the GEdgeUse
};

#endif
