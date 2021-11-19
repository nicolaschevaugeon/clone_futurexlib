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
* $Id: GShell.h,v 1.9 2005/02/02 00:08:55 acbauer Exp $
* Description:

*/
#ifndef H_GShell
#define H_GShell

#include "GEntity.h"

template <class T>
class SSList;
class GRegion;
class GVertex;
class GFace;

/** A Shell. A set of faces uses that define the boundary of a region.
 */
class GShell : public AOMD::GEntity
{
  public:
   GShell(SGModel *m, GRegion *region, const SSList<GFace *> &faces, const SSList<int> &dirs);
   ~GShell() override;

   GeoRep *geometry() override;

   TopoType::Value type() const override;
   int dim() const override;

   AOMD::SBoundingBox3d bounds() const override;

   /** Get the region that this shell is attached to. */
   GRegion *region() const;

   /** Get iterator initialized to first face use of shell. */
   SSListCIter<GFaceUse *> firstFaceUse() const;

   int inClosure(GEntity *ent) const override;

   // these should not be called by anything other than GEntity dervied
   // classes
   GShell(SGModel *m, GRegion *region);  // for reading from file
   GFaceUse *appendFace(GFace *face, int dir);
   void setRegion(GRegion *region);
   void removeFaceUse(GFaceUse *fu);

   void write(ostream &out);
   void read(istream &in);

  protected:
   GRegion *Region;
   SSList<GFaceUse *> FaceUses;
};

#endif
