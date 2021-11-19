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
#ifndef H_GEdgeUsePair
#define H_GEdgeUsePair

#include "GEntity.h"
#include "GEdgeUse.h"

class GLoopUse;
class GEdge;
class GVertex;
class GVertexUse;

/** An edge use. This edge use is a bit different than what is defined for
  a full radial edge representation. Essentially pairs of radial edge uses
  are condensed into one of these edges uses. */
class GEdgeUsePair : public AOMD::GEntity {
public:
  GEdgeUsePair(SGModel *m, GEdge *e, GLoopUse *l, int dir, GVertex *v0, GVertex *v1);
  ~GEdgeUsePair() override;

  // GEntity stuff
  GeoRep * geometry() override;
  AOMD::SBoundingBox3d bounds() const override;
  int inClosure(GEntity *ent) const override; // if ent is in closure of this 
  TopoType::Value type() const override;
  int dim() const override;

  /** Get the owning edge. */
  GEdge *edge() const;

  /** Get the shell. */
  GShell *shell() const;

  /** Get the loop use. */
  GLoopUse *loopUse(int dir) const;

  /** Get the vertex use on the given end. */
  GVertexUse *vertexUse(int dir) const;
  
  /** Return the edge use in the given direction. */
  GEdgeUse *side(int which);

  /** Tell which side the given use side is on. */
  int whichSide(const GEdgeUse *eus) const;
  
  /** Get the other edge use side. */
  GEdgeUse *otherSide(GEdgeUse *eus);

  void write(std::ostream &out);
  void read(std::istream &in);
  //
  GEdgeUsePair(SGModel *m, GEdge *e);
  void setLoopUse(GLoopUse *lu, int dir);
  void setVertexUse(int dir, GVertexUse *vu);
  void mergeWith(GEdgeUsePair *);
  void replaceOwner(GEdge *enew);

  GEdgeUsePair * split(GVertex *v, GEdge *e);
  // debugging
  void info() const;
protected:

private:
  GEdge *OwningEdge;
  GEdgeUse PosSide;
  GEdgeUse NegSide;
};

#endif
