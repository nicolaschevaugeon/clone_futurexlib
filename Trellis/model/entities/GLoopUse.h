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
 /****************************************************************
 * $Id: GLoopUse.h,v 1.13 2005/02/02 00:08:55 acbauer Exp $
 * Description: 

*/
#ifndef H_GLoopUse
#define H_GLoopUse

//#include "SList.h"
#include "GEntity.h"

template<class T> class SSList;
class GEdgeUse;

/** A loop use. Basically a set of edges uses that define a loop on a face. 
  There is no loop that this is a use of at this point. */
class GLoopUse : public AOMD::GEntity {
public:
  GLoopUse(SGModel *m, GFaceUse *, const SSList<GEdge*> &edges, const SSList<int> &dirs);
  ~GLoopUse() override;

  GeoRep * geometry() override;

  TopoType::Value type() const override;
  int dim() const override;

  /** Get the number of edges uses in the loop.  */
  virtual int numEdges() const;

  /** Get iterator initialized to first edge use in loop. */
  SSListCIter<GEdgeUse*> firstEdgeUse() const;

  SSList<GEdge*> edges() const override;

  /** Get face use attached to this loop. */
  virtual GFaceUse* faceUse() const;

  int inClosure(GEntity *ent) const override; // if ent is in closure of this loopuse

  AOMD::SBoundingBox3d bounds() const override;

  // these should not be called by anything other than GEnitity dervied
  // classes
  GLoopUse(SGModel *m, const SSList<GEdgeUse*> & edges, GFaceUse *faceUse, GLoopUse *old);
  GLoopUse(SGModel *m, GFaceUse *); // for reading
  void setMate(GLoopUse *m); // set the loop use mate for this use
  void appendEdge(GEdge *e, int dir);
  void setFaceUse(GFaceUse *f);
  void mergeVertexUses();

  void replace(GEdgeUse *eu, GEdgeUse *newUse);
  void addBefore(GEdgeUse *eu, GEdgeUse *newUse);
  void addAfter(GEdgeUse *eu, GEdgeUse *newUse);
  GLoopUse * split(GEdge *e, int dir, GEdgeUse *after0, GEdgeUse *after1);
  void merge(GEdge *e, int dir, GEdgeUse *after0, GEdgeUse *after1);
  void write(std::ostream &out);
  void read(std::istream &in);

  void info() const;
protected:

  SSList<GEdgeUse *> Edges;
  GFaceUse * FaceUse;
  GLoopUse * d_mate;
};


#endif

