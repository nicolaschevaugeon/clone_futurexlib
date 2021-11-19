/**************************************************************************** 

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written 
   and maintained by the Scientific Computation Research Center (SCOREC) at 
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*****************************************************************************/
#ifndef _MTET_H_
#define _MTET_H_

#include "mRegion.h"
class GEntity;

namespace AOMD {

class mEdge;
class mVertex;
class mFace;
/**
  The mTet class implements specific
  templates for tets
*/
class mTet : public mRegion
{
public:
  ~mTet() override;
  mTet ();
  /// Special constructor with 4 vertices
  mTet (mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, GEntity *classification);
  /// Special constructor with 6 edges
  mTet (mEdge *e1, mEdge *e2, mEdge *e3, mEdge *e4, mEdge *e5, mEdge *e6, GEntity *classification);
  /// Special constructor with 4 faces
  mTet (mFace *f1, mFace *f2, mFace *f3, mFace *f4, GEntity *classification);
  /// Type is mEntity::TET
  inline mEntity::mType getType() const override {return TET;}
  /// template members
  int getNbTemplates (int what) const override;
  /// template members
  mEntity* getTemplate (int ith , int what , int with) const override;
  /// Debug stuff
  void print() const override;
  /// Get the dimension of the entity through a static function
  /// Useful for template algorithms
  inline static int getDim () {return 3;}
};

} // end of namespace

#endif 
