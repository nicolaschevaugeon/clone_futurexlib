/******************************************************************************** 

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

*********************************************************************************/
#include <fstream>
#include <iostream>
#include <cstdio>
#include "mMesh.h"
#include "mEdge.h"
#include "mVertex.h"
#include "mFace.h"
#include "mEntity.h"
#include "AOMD_OwnerManager.h"
#include "ParUtil.h"


namespace AOMD {



  void mMesh::bdryLinkSetup()
  {
    theOwnerManager->clearRemoteCopies();
    int n = size(2)?2:1;
    updateCommonBdryEntities(n);

    for(int dim=0;dim<=n;dim++)
      {
	for(iter it = begin(dim); it != end(dim);++it)
	  {
	    mEntity *e = *it;
	    std::list<mEntity*> mirrors;
	    if(lookupForMirrorEntities (e,mirrors))
	      {
		for( std::list<mEntity*>::const_iterator it2 = mirrors.begin() ;
		     it2 != mirrors.end(); ++it2)
		  {
		    theOwnerManager->addRemoteCopy(e,ParUtil::Instance()->rank(),*it2,true);	  	      
		  }
	      }
	  }
      }  
  }


} // end of namespace

