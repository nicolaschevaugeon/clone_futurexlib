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
  $Id: AttachDataId.h,v 1.2 2005/02/02 00:08:57 acbauer Exp $

File: AttachDataId.h
Description: Header file for attached data management structures
Written by: Jim Teresco
Created: Wed Jun  3 10:10:26 EDT 1998

Modification History:

*/

#ifndef H_AttachDataId
#define H_AttachDataId

#include "SString.h"

class AttachableData;
class AttachDataId;
class AttachDataBase;

typedef void *CallbackClient;

extern "C" {
typedef int (*AttachCallbackProc)(
    AttachableData *    /* entity on which callback occurred */,
    AttachDataId *      /* the data ID which triggered the callback */,
    int                 /* callback event type */,
    CallbackClient *    /* client data - registered by the application */,
    void *              /* call data - specific to each event type */
);
}


class AttachDataId {
  friend class AttachDataBase;
  friend class AttachDataManager;
  friend class AttachableData;
public:
  void setCallback(int, AttachCallbackProc, CallbackClient);
  void removeCallback(int);

  // this should not be public, but will not seem to work otherwise -
  // problem making the templated class AttachedDataV a friend...
  inline void setNeedVirtual(int value) { need_virtual=value; }
protected:
  AttachDataId(const SString &, int);
  virtual ~AttachDataId();
  inline void setID(int idpos) { the_id=idpos; }
  inline int id() { return the_id; }
  virtual void callCallback(int, AttachableData *, void *) = 0;
  int in_callbacks; // this should be a single bit if we have one handy
  int need_virtual; // also only needs one bit
  AttachCallbackProc *functions;
  CallbackClient *client_data;

private:
  int the_id;
  SString idstr;
  int numCBEvent;  // depends on the derived class
};


#endif
