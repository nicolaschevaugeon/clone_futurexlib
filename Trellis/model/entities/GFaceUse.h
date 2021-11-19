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

#ifndef H_GFaceUse
#define H_GFaceUse

#include "GEntity.h"
#include "SSList.h"

/** A face use. Defined as a set of loops. Each face has two of these,
  one on each side. */
class GFaceUse : public AOMD::GEntity {
public:
  GFaceUse(SGModel *m, GFace * face);
  GFaceUse(SGModel *m, GFace * face, GLoopUse * loop, GShell *shell);

  GeoRep * geometry() override;
  TopoType::Value type() const override;
  int dim() const override;

  AOMD::SBoundingBox3d bounds() const override;

  /** Get iterator initialized to first loop use in this face use. */
  SSListCIter<GLoopUse*> firstLoopUse() const;

  /** Get shell use using this face use. */
  GShell * shell() const;
  //SSList<GLoopUse*> loops() const;
  SSList<GEdge *> edges() const override;

  /** Get face assoicated with this face use. */
  GFace * face() const;

  /** Get the direction of this face use relative to the face. */
  int dir() const;

  int inClosure(GEntity *ent) const override; // if ent is in closure of this faceuse

  GLoopUse * insertEdge(GEdge *e);
  void addLoopUse(GLoopUse * lu);
  void setShell(GShell *shell);
  void write(std::ostream &out);  
  void read(std::istream &in);

  // debugging
  void info() const;
private:
  SSList<GLoopUse*> Loops;
  GFace * Face;
  GShell *Shell;
};
#endif
