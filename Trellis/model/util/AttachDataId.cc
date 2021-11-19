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
  $Id: AttachDataId.cc,v 1.2 2005/02/02 00:08:57 acbauer Exp $

File: AttachDataId.cc
Description: attached data management structures
Written by: Jim Teresco
Created: Thu Jun  4 09:35:17 EDT 1998

Modification History:

*/

#include "AttachDataId.h"

AttachDataId::AttachDataId(const SString &id, int numcbevent) {

  in_callbacks=0;
  need_virtual=0;
  the_id=-1;
  idstr=id;
  numCBEvent=numcbevent;
  functions=new AttachCallbackProc[numCBEvent];
  client_data=new CallbackClient[numCBEvent];
  for (int i=0; i<numCBEvent; i++) {
    functions[i]=nullptr;
    client_data[i]=nullptr;
  }
}

AttachDataId::~AttachDataId() {

  delete [] functions;
  delete [] client_data;
}

void AttachDataId::setCallback(int event, AttachCallbackProc proc,
			       CallbackClient client) {

  functions[event]=proc;
  client_data[event]=client;

}

void AttachDataId::removeCallback(int event) {

  functions[event]=nullptr;
  client_data[event]=nullptr;
}

void AttachDataId::callCallback(int event, AttachableData *ent, 
				void *call_data) {
  in_callbacks=1;

  if (functions[event]) {
    (*functions[event])(ent,this,event,&client_data[event],call_data);
  }
  in_callbacks=0;
}


