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

#ifndef _M_FACE_H_
#define _M_FACE_H_
#include "mEntity.h"

  class GEntity;
  
namespace AOMD {


  class mVertex;
  class mEdge;
  /// The Face is implemented generally i.e a face can have N edges
  /// The assumption is that the number
  /// of edges is equal to the number of
  /// vertices and that 2 consecutive edges
  /// share one vertex, that's the template
  /// Special constructors for triangles & quads
  class mFace : public mEntity 
    {
      /// for sppeeding (QUAD or TRI or FACE)
      mType typ;
    public:
      ~mFace() override;
      mFace ();
      /// Special constructor for creating a triangle with vertices
      mFace (mVertex *v1, mVertex *v2, mVertex *v3, GEntity *classification);
      /// Special constructor for creating a triangle with edge
      mFace (mEdge *v1, mEdge *v2, mEdge *v3, GEntity *classification, int *dir = nullptr);
      /// Special constructor for creating a quad with vertices
      mFace (mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, GEntity *classification);
      /// Special constructor for creating a quad with edges
      mFace (mEdge *v1, mEdge *v2, mEdge *v3, mEdge *v4, GEntity *classification, int *dir = nullptr);
      /// getLevel returns 2
      int getLevel()const override;
      /// getType could return QUAD, TRI or FACE
      inline mEntity::mType getType ()const override{return typ;}
      /// Get the common vertex of 3 faces
      mVertex *commonVertex (mFace *f1, mFace *f2);
      /// Get the common edge of 2 faces
      mEdge *commonEdge (mFace *f1);
      /// Template members
      int getNbTemplates (int what) const override;
      /// Template members
      mEntity* getTemplate (int ith , int what , int with) const override;
      /// Debug stuff
      void print() const override;
      /// Get the dimension of the entity through a static function
      /// Useful for template algorithms
      inline static int getDim () {return 2;}
  
    };

} // end of namespace

#endif 

