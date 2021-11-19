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
 * $Id: GRegion.h,v 1.10 2005/02/02 00:08:55 acbauer Exp $
 * Description: 

*/
#ifndef H_GRegion
#define H_GRegion

//#include "SList.h"
#include "GEntity.h"
//#include "GFace.h"

class GShell;
class GFace;
template<class T> class SSList;

/** A model region. */
class GRegion : public AOMD::GEntity {
public:
  GRegion(SGModel *model, int tag, const SSList<GFace*> & faces, const SSList<int> &dirs);

  ~GRegion() override;

  TopoType::Value type() const override;
  int dim() const override;

  SSList<GRegion*> regions() const override;
  SSList<GFace*> faces() const override;
  SSList<GEdge*> edges() const override;
  SSList<GVertex*> vertices() const override;

  /** Get iterator initialized for first shell in the region definition. */
  SSListCIter<GShell*> firstShell() const;

  int inClosure(GEntity *ent) const override; 
  int classifyPoint(const AOMD::SPoint3 &pt, int chkb, GEntity **bent) const override;
  AOMD::SBoundingBox3d bounds() const override;

  // these should not be called by anything other than GEntity dervied
  // classes
  GRegion(SGModel *model, int tag);
  void appendShell(GShell *shell);
  void addShell(const SSList<GFace*> & faces, const SSList<int> &dirs);
  void write(std::ostream &out);
  void read(std::istream &in);

protected:
  SSList<GShell*> Shells;

};

#endif

