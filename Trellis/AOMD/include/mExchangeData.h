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
#ifndef _AOMD_EXCHANGEDATA_H_
#define _AOMD_EXCHANGEDATA_H_

#include "mMesh.h"
#include "AOMD_SharedInfo.h"

namespace AOMD {

class mEntity;

/**
  AOMD_DataExchanger is a class that allow user to exchange
  data's between partitions and/or periodic boundaries. User
  has to provide a tag for non synchronous communications. 
  Basically, choose a tag that is higher that 200 for being
  sure that it does not interact with internal messages.
  User has to allocate buffers using AP_alloc from autopack.
  this should be changed... We already did the change in 
  load balancing callbacks.
  User will recieve its datas. If user wants to know the
  counterpart of mesh entity *e in the other partition,
  he has to sent the remote pointer si.getRemotePointer().
  This also should be automatized.
*/
  class AOMD_DataExchanger
  {
  public :
    virtual  ~AOMD_DataExchanger() = default;;
  /// get a tag
    virtual int tag() const = 0;
    /// user sends a message related to mesh entity *e to proc si.pid() to the counterpart si.getRemotePointer().
    virtual void * AP_alloc_and_fill_buffer (mEntity *e, AOMD_SharedInfo &si, int tag) = 0;
    /// user recieve data *buf form proc pid.
    virtual void receiveData (int pid, void *buf) = 0; 
  };
} // end of namespace 

#endif
