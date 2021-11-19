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
#ifndef _MREGION_H_
#define _MREGION_H_

#include "mEntity.h"

namespace AOMD {

class mEdge;
class mFace;
/**
   mRegion is the base class for all 3D elements
   It's impossible to find a common template for
   all 3D elements like we did with faces. The only
   thing that is common is that the dimension is 3.
   
   Some members are also common to all regions.
*/
  class mRegion  : public mEntity
  {
  public:
    mRegion(): mEntity(){}
    ~mRegion() override= default;
    /// Dimension is 3
    inline int getLevel()const override{return 3;};
    /// Gives the common edge to 3 regions
    mEdge *commonEdge (mRegion *r1, mRegion *r2);
    /// Gives the common face to 2 regions
    mFace *commonFace (mRegion *r);
    /// Get the dimension of the entity through a static function
    /// Useful for template algorithms
    inline static int getDim () {return 3;}
  }; 
} // end of namespace

#endif 
