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
  $Id: ModelDataId.h,v 1.8 2005/02/02 00:08:55 acbauer Exp $

File: ModelDataId.h
Description: Header file for attached data management structures
Written by: Jim Teresco
Created: Fri Jun 19 15:49:40 EDT 1998 

Modification History:

*/

#ifndef H_ModelDataId
#define H_ModelDataId

#include "AttachDataId.h"
//#include "utilfwd.h"


class SModelMember;
class ModelDataId;

typedef void *CallbackClient;

extern "C" {
typedef int (*ModelCallbackProc)(
    SModelMember *  /* entity on which callback occurred */,
    int             /* callback event type */,
    CallbackClient  /* client data - registered by the application */,
    void *          /* call data - specific to each event type */,
    ModelDataId *    /* the data ID which triggered the callback */
);
}

class ModelDataId: public AttachDataId {
  friend class AttachDataBase;
  friend class ModelDataManager;
  friend class SModelMember;
public:
  enum ModelCBEvent {
    Create=0,
    Delete,
    Modify,
    NumModelCBEvent
  };

protected:
  ModelDataId(const SString &);
  ~ModelDataId() override;
  void callCallback(int, AttachableData *, void *) override;
};


#endif
