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
 /***************************************************************
 * $Id: TRegion.h,v 1.9 2005/02/02 00:08:57 acbauer Exp $
 * Description: 
 * Modification 9/15/98 - added the tagId member data that is returned in the
 * case of that RealRegion = Null Pointer 

*/

#ifndef H_TRegion
#define H_TRegion

//#include "SList.h"
#include "GRegion.h"

template<class T> class SSList;
class GFace;
class GeoRep;

/** A topological region. May represent a copy of a region in another model. */
class TRegion : public GRegion {
public:
  TRegion(SGModel *model, int tag, GRegion *r,const SSList<GFace*> & faces, 
	  const SSList<int> &dirs);
  TRegion(SGModel *model, int tag);
  ~TRegion() override;

  int containsPoint(const AOMD::SPoint3 &pt) const override;

  GeoRep * geometry() override;

  double tolerance() const override;

protected:
  GRegion *RealRegion;
};


#endif
