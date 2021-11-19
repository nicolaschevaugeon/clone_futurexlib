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
#include "AttachableData.h"
#include "TList.h"
#include "AttachData.h"
#include "AttachDataManager.h"
#include "AttachDataId.h"

AttachableData::AttachableData()
  : Data(nullptr), Data1(nullptr)
{}

AttachableData::~AttachableData()
{
  if(Data)
    TList_delete(Data);

  // the section added by jie starts
  AttachDataBase *next;
  for (AttachDataBase *current=Data1; current; current=next) {
    next =current->next;
    delete current;
  }
  // the section added by jie ends
}

void AttachableData::setDataP(const char *tag, void *d)
{
  int ntag;
  ntag = (((int)tag[0])<<24) | (((int)tag[1])<<16) | (((int)tag[2])<<8) | tag[3];
  if( !Data)
    Data = TList_new();
  TList_append(Data,ntag,d);

}

void AttachableData::setDataI(char *tag, long int d)
{
  int ntag;

  ntag = (((int)tag[0])<<24) | (((int)tag[1])<<16) | (((int)tag[2])<<8) | tag[3];
  if( !Data )
    Data = TList_new();

  TList_append(Data,ntag,(void*)d);

}

void * AttachableData::dataP(const char*tag) const
{
  TListNode *node;
  int ntag;

  ntag = (((int)tag[0])<<24) | (((int)tag[1])<<16) | (((int)tag[2])<<8) | tag[3];
  if(Data){
    node = TList_itemInList(Data,ntag);
    if(node)
      return node->data;
  }
  return nullptr;

}

long int AttachableData::dataI(char *tag) const
{
  TListNode *node;
  int ntag;

  ntag = (((int)tag[0])<<24) | (((int)tag[1])<<16) | (((int)tag[2])<<8) | tag[3];
  if(Data){
    node = TList_itemInList(Data,ntag);
    if(node)
      return (long int)(node->data);
  }
  return 0;

}

void AttachableData::removeData(char *tag)
{ 
  int ntag;

  ntag = (((int)tag[0])<<24) | (((int)tag[1])<<16) | (((int)tag[2])<<8) | tag[3];

  if(Data)
    TList_remItem(Data, ntag);

}

int AttachableData::modifyDataP(char *tag, void * data)
{
  TListNode *node;
  int ntag;

  if(Data){
    ntag = (((int)tag[0])<<24) | (((int)tag[1])<<16) | (((int)tag[2])<<8) | tag[3];
    node = TList_itemInList(Data, ntag);
    if(node){
      node->data = data;
      return 1;
    }
  }
  return 0;
}

int AttachableData::modifyDataI(char *tag, long int data)
{
  TListNode *node;
  int ntag;

  if(Data){
    ntag = *(int *)tag;
    node = TList_itemInList(Data, ntag);
    if(node){
      node->data = (void *)data;
      return 1;
    }
  }
  return 0;
}

// the section added by jie starts
void AttachableData::attachData(AttachDataBase *attach) {

  if (!Data1) {
    Data1=attach;
    attach->next=nullptr;
    return;
  }

  // insert in the embedded list, in ascending id order
  if (Data1->id() > attach->id()) {
    // it goes first
    attach->next=Data1;
    Data1=attach;
    return;
  }
  AttachDataBase *trail, *current;
  
  trail=Data1;
  current=Data1->next;
  while (current) {
    if (current->id() > attach->id()) {
      attach->next=current;
      trail->next=attach;
      return;
    }
    trail=current;
    current=current->next;
  }
  
  // it goes at the end
  trail->next=attach;
  attach->next=nullptr;
}

void AttachableData::deleteData(AttachDataId *id, int warn_exist)
{
  delete detachData(id, warn_exist);
}

AttachDataBase * AttachableData::detachData(AttachDataId *id, int warn_exist) 
{
  if (!Data1) {
    return nullptr;
  }
  AttachDataBase *ret;

  if (Data1->id()==id->id()) {
    // remove first entry in list
    ret = Data1;
    Data1=Data1->next;
    ret->next=nullptr;
    return ret;
  }
  AttachDataBase *trail=Data1;
  AttachDataBase *current=Data1->next;
  while (current) {
    if (current->id()==id->id()) {
      ret = current;
      trail->next=current->next;
      current->next=nullptr;
      return ret;
    }
    trail=current;
    current=current->next;
  }

  return nullptr;
}

AttachDataBase *AttachableData::getData(AttachDataId *id) const {

  for (AttachDataBase *current=Data1; current; current=current->next) {
    if (current->id()==id->id())
      return current;
  }
  return (AttachDataBase *)nullptr;
}

void AttachableData::callCallbacks(int event, void *call_data) {
  // create is a special case since no data has yet been attached
  // so we call callbacks for all known tags to give applications
  // the chance to attach anything they want
  if (event == 3) {  // this assumes create is always 3 - this is bad.
    AttachDataIdIter *iter=dataManager()->firstAttachDataId();
    AttachDataId *dataid;
    while ((dataid=iter->next())) {
      // call create callbacks associated with this one
      dataid->callCallback(event, this, call_data);
    }
    delete iter;
  }
  else {
    // loop over data attached to this entity and call the callbacks
    for (AttachDataBase *current=Data1; current; current=current->next) {
      AttachDataId *dataid=current->attachdataid();
      dataid->callCallback(event, this, call_data);
    }
  }
}

void AttachableData::attachDataInt(AttachDataId *id, int value)
{ attachData(new AttachedData<int>(id,value)); }

void AttachableData::attachDataDbl(AttachDataId *id, double value)
{ attachData(new AttachedData<double>(id,value)); }

void AttachableData::attachDataPtr(AttachDataId *id, void * value)
{ attachData(new AttachedData<void *>(id,value)); }

int AttachableData::getDataInt(AttachDataId *id, int *value) const
{
  AttachedData<int> *attached=(AttachedData<int> *)getData(id);
  *value = attached ? attached->data : 0;
  return attached!=nullptr;
}

int AttachableData::getDataDbl(AttachDataId *id, double *value) const
{
  AttachedData<double> *attached=(AttachedData<double> *)getData(id);
  *value = attached ? attached->data : 0.0;
  return attached!=nullptr;
}

int AttachableData::getDataPtr(AttachDataId *id, void **value) const
{
  AttachedData<void *> *attached=(AttachedData<void *> *)getData(id);
  *value = attached ? attached->data : nullptr;
  return attached!=nullptr;
}

void AttachableData::modifyDataInt(AttachDataId *id, int value)
{
  AttachedData<int> *attached=(AttachedData<int> *)getData(id);
  if(attached)
    attached->data = value;
  else
    attachDataInt(id,value);
}

void AttachableData::modifyDataDbl(AttachDataId *id, double value)
{
  AttachedData<double> *attached=(AttachedData<double> *)getData(id);
  if(attached)
    attached->data = value;
  else
    attachDataDbl(id,value);
}

void AttachableData::modifyDataPtr(AttachDataId *id, void * value)
{
  AttachedData<void *> *attached=(AttachedData<void *> *)getData(id);
  if(attached)
    attached->data = value;
  else
    attachDataPtr(id,value);
}
// the section added by jie ends
