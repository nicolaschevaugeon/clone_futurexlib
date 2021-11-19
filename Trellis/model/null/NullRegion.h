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
#ifndef H_NullRegion
#define H_NullRegion

 /***************************************************************
 * $Id: NullRegion.h,v 1.8 2005/02/02 00:08:56 acbauer Exp $
 * Description: 

*/
//#include "SList.h"
#include "GRegion.h"
#include "NullEntity.h"

template<class T> class SSList;
class GFace;
class NullFace;
class GeoRep;

/** A "fake" model region. */
class NullRegion : public GRegion, public NullEntity {
public:
  NullRegion(NullModel *m, int tag, const SSList<GFace*> &faces, const SSList<int> &dirs);
  NullRegion(NullModel *m, int tag);
  ~NullRegion() override;

  GeoRep * geometry() override;

  double tolerance() const override;

  int containsPoint(const AOMD::SPoint3 &pt) const override;

protected:

};


#endif
