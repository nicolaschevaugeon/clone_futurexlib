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
#ifndef _MPYRAMID_H_
#define _MPYRAMID_H_
#include "mRegion.h"
class mEdge;
class mVertex;
class mFace;
/**
  mPyramid implement specific templates for the pyramidal
  element.
 */
class mPyramid : public mRegion
{
 public:
  virtual ~mPyramid();
  mPyramid ();
  /// Constructor using 5 vertices
  mPyramid (mVertex *v1, 
	mVertex *v2, 
	mVertex *v3, 
	mVertex *v4, 
	mVertex *v5, 
	GEntity *classification);
  /// Constructor using 8 edges
  mPyramid (mEdge *e1, mEdge *e2, mEdge *e3, mEdge *e4, mEdge *e5, mEdge *e6, 
	mEdge *e7, mEdge *e8, GEntity *classification);
  /// Constructor using 5 face
  mPyramid (mFace *f1, mFace *f2, mFace *f3, mFace *f4, mFace *f5, 
	GEntity *classification);
  /// This is an mEntity::PYRAMID
  inline mEntity::mType getType() const {return PYRAMID;}
  /// Template members
  virtual int getNbTemplates (int what) const;
  /// Template members
  virtual mEntity* getTemplate (int ith , int what , int with) const;
  /// Debug stuff
  virtual void print() const;
};

#endif 
