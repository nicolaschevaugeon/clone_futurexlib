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
  $Id: AttachDataManager.h,v 1.2 2005/02/02 00:08:57 acbauer Exp $

File: AttachDataManager.h
Description: Header file for attached data management class
Written by: Jim Teresco
Created: Wed Jun  3 10:10:07 EDT 1998

Modification History:

*/

#ifndef H_AttachDataManager
#define H_AttachDataManager

#include "AttachDataId.h"
#include "SString.h"

template<class T> class SBlock;
class AttachDataIdIter;

class AttachDataManager {
  friend class AttachDataIdIter;
  friend class AttachableData;
public:
  AttachDataManager();
  virtual ~AttachDataManager();
  void unregisterId(AttachDataId *);
  AttachDataId *lookupTagId(const char *tag);

protected:
  AttachDataIdIter *firstAttachDataId();
  // these members must be called from all processes in the same order
  // they should be called from the public members newMeshDataId, etc
  // in derived classes
  virtual AttachDataId *newAttachDataId(const SString &) = 0;
  void registerAttachDataId(AttachDataId *);

  int register_count;
  SBlock<AttachDataId *> *idlist;
  int initialized;
};

class AttachDataIdIter {
  friend class AttachDataManager;
  friend class AttachableData;
protected:
  AttachDataIdIter(AttachDataManager *man) { current_pos=0; idlist=man->idlist; }
  ~AttachDataIdIter() = default;
  AttachDataId *next();

private:
  SBlock<AttachDataId *> *idlist;
  int current_pos;
};

#endif
