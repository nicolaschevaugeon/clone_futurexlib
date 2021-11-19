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
#ifndef H_GVertex
#define H_GVertex

#include "GEntity.h"
#include "SPoint3.h"

template<class T> class SSList;
class GEdge;
class GEdgeUsePair;
class GVertexUse;
class GVPoint;

/** A model vertex. */
class GVertex : public AOMD::GEntity {
public:
  GVertex(SGModel *m, int tag);
  ~GVertex() override = default;

  TopoType::Value type() const override;
  int dim() const override;

  int inClosure(GEntity *ent) const override; // if ent is in closure of this vertex

  SSList<GRegion*> regions() const override;
  SSList<GFace*> faces() const override;
  SSList<GEdge*> edges() const override;
  SSList<GVertex*> vertices() const override;

  /** Number of vertex uses. */
  int numUses() const;

  /** Get the vertex uses at this vertex. */
  SSList<GVertexUse*> uses() const;  

  /** Get the number of edges using this vertex. If the same edges uses
     a vertex twice (at each end), it is counted twice. */
  int numEdges() const;

  /** Get the classified location of this vertex. */
  virtual GVPoint point() const = 0;

  int classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const override;
  int containsPoint(const AOMD::SPoint3 &pt) const override;

  AOMD::SBoundingBox3d bounds() const override;

  // debugging
  void info() const;

  // these should not be called by anything other than GEnitity dervied
  // classes
  GVertexUse * addEdge(GEdgeUsePair *e);
  void removeUse(GVertexUse *vu);

  void write(std::ostream &out);
  void read(std::istream &in);
  GVertexUse * use(int id);

protected:
  SSList<GVertexUse *> Uses;

};

#endif

