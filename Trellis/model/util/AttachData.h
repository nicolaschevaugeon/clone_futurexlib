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
  $Id: AttachData.h,v 1.3 2005/02/02 00:08:57 acbauer Exp $

File: AttachData.h
Description: Header file for structures to store attached data
Written By: Jim Teresco
Created: Thu Jun  4 11:02:26 EDT 1998

Modification History:

*/

#ifndef H_AttachData
#define H_AttachData

#include "AttachDataId.h"
#include <iostream>

class AttachDataId;

// base class for all data that can be attached
class AttachDataBase {
  friend class AttachableData;
public:
  AttachDataBase(AttachDataId *id) { myid=id; }
  ~AttachDataBase() = default;

protected:
  int id() { return myid->id(); }
  AttachDataId *attachdataid() { return myid; }
  // embedded list pointer
  AttachDataBase *next;
private:
  AttachDataId *myid;
};

// templated class for simple data types (those that have 
// trivial destructors). Since the destructor of AttachDataBase
// isn't virtual, the derived destructor won't get called. This
// is by design to avoid the overhead of have a vtable in the
// base class.
// Really, this needs to be looked at, since the overhead of
// the system new operator may be more than what we would
// pay for a vtable pointer and without it we can't use our
// zero overhead allocator
template<class T>
class AttachedData: public AttachDataBase {
public:
  AttachedData(AttachDataId *id, T &value):AttachDataBase(id) { data=value; }
  AttachedData(AttachDataId *id):AttachDataBase(id) {}
  ~AttachedData() = default;
  T data;
  //void setData(T value) { current_val=value; }
  //T &data() { return current_val; }
};

// base class for attached data that have a non trivial destructor
// this needs to be fixed, should have a class AttachedDataVirtual
// derived from AttachData, then a templated class derived from that.
template<class T>
class AttachedDataV: public AttachedData<T> {
public:
  AttachedDataV(AttachDataId *id, T &value);
  AttachedDataV(AttachDataId *id);
  virtual ~AttachedDataV();
};

template<class T>
inline AttachedDataV<T>::AttachedDataV(AttachDataId *id, T &value):AttachedData<T>(id,value)  {

  id->setNeedVirtual(1);
}

template<class T>
inline AttachedDataV<T>::AttachedDataV(AttachDataId *id):AttachedData<T>(id) {

  id->setNeedVirtual(1);
}

template<class T>
AttachedDataV<T>::~AttachedDataV() {

  // just for now - debugging
  std::cerr << "Called the virtual destructor AttachedDataV<T>::~AttachedDataV()\n";
}

#endif
