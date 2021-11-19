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
#ifndef  _AOMD_EXCHANGEDATAGENERIC_H_
#define  _AOMD_EXCHANGEDATAGENERIC_H_

#include "mExchangeData.h"
#include "AOMD_OwnerManager.h"
#include "mMesh.h"
#include "ParUtil.h"


namespace AOMD {
  

  template <class Iterator> 
  void exchangeDataGeneric ( const Iterator &beg, const Iterator &end , AOMD_DataExchanger &de )
  {
    Iterator it = beg;
    
    for( ; it != end ; ++it)
      {	      
	mEntity *e = (*it).first;
	AOMD_SharedInfo si = (*it).second;    
	void *buf = de.AP_alloc_and_fill_buffer (e, si, de.tag());
	if(buf)
	  {
	    de.receiveData (ParUtil::Instance()->rank(), buf);
	    free(buf);
	  }
      }   
  }

  

} // end of namespace
#endif // end define

