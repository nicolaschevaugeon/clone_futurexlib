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
  $Id: AttachDataManager.cc,v 1.2 2005/02/02 00:08:57 acbauer Exp $

File: AttachDataManager.cc
Description: attached data management class
Written by: Jim Teresco
Created: Wed Jun  3 11:06:12 EDT 1998

Modification History:

*/

#include "AttachDataManager.h"
#include "SBlock.h"

#define ALLOCSIZE 5


AttachDataManager::AttachDataManager() {

  int i;
  if (!initialized) {
    initialized=1;
  }
  else {
    // initialized twice - error 
  }
  register_count=0;
  idlist=new SBlock<AttachDataId *>(ALLOCSIZE);
  for (i=0; i<idlist->size(); i++) 
    (*idlist)[i]=nullptr;
}

AttachDataManager::~AttachDataManager() {

  // derived classes need to delete the AttachDataId's
  // delete any AttachDataId's still registered
  //for (int i=0; i<idlist->size(); i++) {
  //  if ((*idlist)[i]) delete (*idlist)[i];
  //}
  delete idlist;
}

void AttachDataManager::registerAttachDataId(AttachDataId *attach) {
  int newpos=-1;

  if (register_count == idlist->size()) {
    // all slots full, make more slots
    idlist->size(idlist->size()+ALLOCSIZE);
    for (int i=register_count; i<idlist->size(); i++)
      (*idlist)[i]=nullptr;
    newpos=register_count;
  }
  else {
    // space available - find a free slot
    for (int i=0; i<idlist->size(); i++) {
      if (!(*idlist)[i]) {
	newpos=i;
	break;
      }
    }
    if (newpos==-1) {
      //InternalError("logic error - new identifier not found");
    }
  }
  (*idlist)[newpos]=attach;
  attach->setID(newpos);
  register_count++;
}

struct tagid {
  int tag;
  AttachDataId *dataid;
};

AttachDataId *AttachDataManager::lookupTagId(const char *tag) {
  // translation of tags to AttachDataId's is done entirely here
  static SBlock<tagid> tagidlist;

  int ntag=(((int)tag[0])<<24) | (((int)tag[1])<<16) | (((int)tag[2])<<8) | tag[3];
  if (ntag==0) {
    //Error("AttachDataManager::lookupTagId: reserved tag passed in");
  }

  int i;
  // see if it's already there
  for (i=0; i<tagidlist.size(); i++) {
    if (ntag == tagidlist[i].tag) return tagidlist[i].dataid;
  }
 
  // not in the list, add it and return it
  int newpos=-1;
  for (i=0; i<tagidlist.size(); i++) {
    if (tagidlist[i].tag == 0) {
      newpos=i;
      break;
    }
  }
  if (newpos == -1) {
    newpos=tagidlist.size();
    tagidlist.size(tagidlist.size()+5);
    for (i=newpos; i<tagidlist.size(); i++) {
      tagidlist[i].tag=0;
      tagidlist[i].dataid=nullptr;
    }
  }
  tagidlist[newpos].tag=ntag;
  tagidlist[newpos].dataid=newAttachDataId(SString(tag));
  return tagidlist[newpos].dataid;
}

void AttachDataManager::unregisterId(AttachDataId *id) {

  for (int i=0; i<idlist->size(); i++) {
    if ((*idlist)[i] == id) {
      //delete (*idlist)[i];  // derived classes are repsonsible for this
      (*idlist)[i]=nullptr;
      register_count--;
      return;
    }
  }
  //Error("AttachDataId not found");
}

AttachDataIdIter *AttachDataManager::firstAttachDataId() {

  return new AttachDataIdIter(this);
}

AttachDataId *AttachDataIdIter::next() {

  while (current_pos<idlist->size()) {
    if ((*idlist)[current_pos]) {
      current_pos++;
      return (*idlist)[current_pos-1];
    }
    current_pos++;
  }
  return nullptr;
}


