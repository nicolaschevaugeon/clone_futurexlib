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

#include <iostream>
#include "mCommonBdry.h"

namespace AOMD {

  bool mCommonBdryLessThanKey::operator() (mCommonBdry* cb1, mCommonBdry* cb2) const
  {
    if(cb1->getId() < cb2->getId())return true;
    if(cb1->getId() > cb2->getId())return false;
    if(cb1->getLevel() < cb2->getLevel())return true;
    if(cb1->getLevel() > cb2->getLevel())return false;

    return false;
  }

} // end of namespace
