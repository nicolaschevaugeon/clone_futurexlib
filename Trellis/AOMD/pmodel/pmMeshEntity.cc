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

/*<i> ****************************************************************************
 *
 *  AOMD/pmodel/pmMeshEntity.cc
 *  Created by Seegyoung Seol, on Mon Dec 08 2003, 10:26:37 EDT
 *
 *  File Content: mEntity operators dealing with PClassification
 *
 *************************************************************************** </i>*/


#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <cassert>
#include <vector>

#include "mEntity.h"
#include "pmEntity.h"
#include "ParUtil.h"

using std::copy;
using std::map;
using std::cout;
using std::endl;
using std::set;
using std::vector;

namespace AOMD {

int mEntity::getOwner()
{
  return 0;

}


void mEntity::setPClassification(pmEntity* pe)
{

}

pmEntity* mEntity::getPClassification()
{

  return (pmEntity*)nullptr;

}



} /* end of namespace */ 
