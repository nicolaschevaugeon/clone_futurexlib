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

#ifndef _M_EDGE_H_
#define _M_EDGE_H_
#include "mEntity.h"

  class GEntity;
  
  namespace AOMD {


  class mVertex;
  /**
     The mEdge Class is the base class for
     edges i.e. entities of dimension 1.

     Edges have always 2 vertices, that's
     the only template that has to be provided.  
  */
  class mEdge : public mEntity 
    {
    public:
      ~mEdge() override;
      /// Special constructor taking 2 vertices as input
      mEdge (mVertex *v1, mVertex *v2,GEntity *classification);
      /// getLevel returns 1 of course
      int getLevel()const override;
      /// getType returns mEntity::EDGE
      inline mEntity::mType getType() const override{return EDGE;}
      /// a shortcut for getting vertices easily
      inline mVertex *vertex(int) const;
      /// gives the common vertex of 2 edges
      inline mVertex *commonVertex (mEdge *) const;
      /// templates for an edge are basically none
      /// but function has to be overwritten
      int getNbTemplates (int what) const override;
      /// debug stuff
      void print() const override;
      /// Get the dimension of the entity through a static function
      /// Useful for template algorithms
      inline static int getDim () {return 1;}
      void setVertices (mVertex *v1, mVertex *v2);
    };

  inline mVertex *mEdge::vertex(int i) const
    {
      return (mVertex*)theAdjacencies[0]->get(i);
    }

  inline mVertex *mEdge::commonVertex(mEdge *other) const
    {
      if(other->vertex(0) == vertex(0))return vertex(0);
      if(other->vertex(1) == vertex(0))return vertex(0);
      if(other->vertex(0) == vertex(1))return vertex(1);
      if(other->vertex(1) == vertex(1))return vertex(1);
      return nullptr;
    }
} // end of namespace

#endif 

