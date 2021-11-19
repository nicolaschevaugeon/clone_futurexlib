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
#ifndef _MPRISM_H_
#define _MPRISM_H_
#include "mRegion.h"
namespace AOMD {

  class mEdge;
  class mVertex;
  class mFace;
  /**
     mPrism implement specific templates for the prismatic (or wedge)
     element.
  */
  class mPrism : public mRegion
    {
    public:
      ~mPrism() override;
      mPrism ();
      /// Constructor using 6 vertices
      mPrism (mVertex *v1, 
	      mVertex *v2, 
	      mVertex *v3, 
	      mVertex *v4, 
	      mVertex *v5, 
	      mVertex *v6,
	      GEntity *classification);
      /// Constructor using 9 edges
      mPrism (mEdge *e1, mEdge *e2, mEdge *e3, mEdge *e4, mEdge *e5, mEdge *e6, 
	      mEdge *e7, mEdge *e8, mEdge *e9, GEntity *classification);
      /// Constructor using 5 face
      mPrism (mFace *f1, mFace *f2, mFace *f3, mFace *f4, mFace *f5, 
	      GEntity *classification);
      /// This is an mEntity::PRISM
      inline mEntity::mType getType() const override {return PRISM;}
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
}

#endif 
