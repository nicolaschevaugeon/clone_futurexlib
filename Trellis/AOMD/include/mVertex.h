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


#ifndef _M_VERTEX_H_
#define _M_VERTEX_H_

#include "mPoint.h"
#include "mEntity.h"

  class GEntity;
  
  namespace AOMD {

  class mIdGenerator;
  /**
     mVertex is the vertex class i.e. entities of
     dimsnsion 0.
  */

  class mVertex : public mEntity 
    {
  
    protected:
      /// A point (position)
      Trellis_Util::mPoint p;
      /// a random value computed from vertex iD
      /// we store this value because random values are
      /// expensive to compute
      int RAND;
    public:      
      /// Constructor taking an id as input
      mVertex(int theId, const Trellis_Util::mPoint & ,GEntity *classif);
      /// Constructor taking an idGenerator as input
      mVertex(mIdGenerator &theIdGenerator , const Trellis_Util::mPoint &, GEntity *classif);
      /// returns a copy of the point
      inline Trellis_Util::mPoint point() const {return p;};
      /// only for meshsim compatibility (DO NOT USE THAT)
      inline Trellis_Util::mPoint* ppoint() {return &p;};
      /// returns RAND
      int getRAND();
      /// delete the id in the id generator (destruction)
      void deleteId(mIdGenerator &theIdGenerator);
      ~mVertex() override;
      /// Dimsnsion is 0
      int getLevel() const override;
      /// Type is mEntity::VERTEX
      mEntity::mType getType() const override{return VERTEX;}
      /// Debug stuff
      void print() const override;
      /// Get the dimension of the entity through a static function
      /// Useful for template algorithms
      inline static int getDim () {return 0;};
      /// move the vertex to another location
      inline void move (const Trellis_Util::mPoint &pt) { p = pt;}
      /// change the iD of the Vertex (dangerous !!!)
      inline void setId( int id ) {iD = id;srand(iD);RAND = abs(rand()) % 100000000;}
      inline void setIdNoRand( int id ) {iD = id;}
    };

  class VertexLexicographicLessThan 
    {
      double EPS;
    public:
      VertexLexicographicLessThan (double eps = 1.e-6); 
      bool operator()(mVertex* ent1, mVertex* ent2) const;
    };

} // end of namespace

#endif 

