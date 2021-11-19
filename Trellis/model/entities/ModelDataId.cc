/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/

/*
  $Id: ModelDataId.cc,v 1.7 2005/02/02 00:08:55 acbauer Exp $

File: ModelDataId.cc
Description: attached data management structures
Written by: Jim Teresco
Created: Thu Jun 18 16:15:54 EDT 1998

Modification History:

*/

#include "ModelDataId.h"
#include "SModelMember.h"
//#include "MessageOut.h"


ModelDataId::ModelDataId(const SString &id) : AttachDataId(id, NumModelCBEvent) {
}

ModelDataId::~ModelDataId() {

  //AttachDataId::~AttachDataId();
}

void ModelDataId::callCallback(int event, AttachableData *ent, 
			       void *call_data) {
//   in_callbacks=1;
//   if (functions[event]) {
//     ModelCallbackProc callback=(ModelCallbackProc)functions[event];
//     (*callback)((SModelMember *)ent,event,client_data[event],call_data,this);
//   }
//   in_callbacks=0;
}


