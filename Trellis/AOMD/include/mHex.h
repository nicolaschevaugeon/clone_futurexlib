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

#ifndef _MHEX_H_
#define _MHEX_H_
#include "mRegion.h"

class GEntity;

namespace AOMD {

class mEdge;
class mVertex;
class mFace;
/**
  mHex implement specific templates for the hexaedral
  element.
 */
class mHex : public mRegion
{
 public:
  ~mHex() override;
  mHex ();
  /// Constructor using 8 vertices
  mHex (mVertex *v1, 
	mVertex *v2, 
	mVertex *v3, 
	mVertex *v4, 
	mVertex *v5, 
	mVertex *v6, 
	mVertex *v7, 
	mVertex *v8, 
	GEntity *classification);
  /// Constructor using 12 edges
  mHex (mEdge *e1, mEdge *e2, mEdge *e3, mEdge *e4, mEdge *e5, mEdge *e6, 
	mEdge *e7, mEdge *e8, mEdge *e9, mEdge *e10, mEdge *e11, mEdge *e12, 
	GEntity *classification);
  /// Constructor using 6 face
  mHex (mFace *f1, mFace *f2, mFace *f3, mFace *f4, mFace *f5, 
	mFace *f6, GEntity *classification);
  /// This is an mEntity::HEX
  inline mEntity::mType getType() const override {return HEX;}
  /// Template members
  int getNbTemplates (int what) const override;
  /// Template members
  mEntity* getTemplate (int ith , int what , int with) const override;
  /// Debug stuff
  void print() const override;
      /// Get the dimension of the entity through a static function
      /// Useful for template algorithms
      inline static int getDim () {return 3;}
};

} // end of namespace

#endif 
