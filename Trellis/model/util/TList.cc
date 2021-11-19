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
#include "TList.h"
#include <cstdlib>

TListNode * TListNode_new()
{ 
  TListNode *object;
  object = new TListNode;
  object->next = nullptr; 
  object->data = nullptr;
  return object;
}

TListNode *TListNode_newItem(int  tag, void *data)
{
  TListNode *object;
  object = new TListNode;
  object->next =nullptr;
  object->tag = tag;
  object->data = data;
  return object;
}

TListNode *TListNode_newCopy( TListNode *nd)
{
  TListNode *object;
  object = new TListNode;
  object->next = nullptr;
  object->tag = nd->tag;
  object->data = nd->data;
  return object;
}

void TListNode_delete(TListNode *object) 
{
  delete object;
}

TListNode * TListNode_assign(TListNode *object,TListNode *node) 
    { object->next = node->next; object->tag=node->tag; return object; }

TListNode * TListNode_next(TListNode *object)
      { return object->next; }

TList * TList_new()
{ 
  TList *object;
  object = new TList;
  object->head = nullptr;
  object->tail = nullptr;
  object->size = 0;
  return object;
} 

  
TList * TList_newCopy(TList *il)
{
  TList *object;
  TListNode *cur;
  TListNode *nd = il->head;

  object = new TList;

  cur = object->head = TListNode_newCopy(nd);
  object->size = 1;
  while(nd->next){
    nd = nd->next;
    cur->next = TListNode_newCopy(nd);
    cur = cur->next;
    object->size ++;
  }
  return object;
}
 
void TList_delete(TList *object)
{
  TListNode *ndNext;
  TListNode *ndCur = object->head;
  while(ndCur){
    ndNext = ndCur->next;
    TListNode_delete(ndCur);
    ndCur = ndNext;
  }
  free(object);
}


/*TList & TList::operator=(TList & il)*/
TList * TList_copy(TList *object,TList *il)
{
  TListNode *cur;
  TListNode *nd = il->head;
  cur = object->head = TListNode_newCopy(nd);
  object->size = 1;
  while(nd->next){
    nd = nd->next;
    cur->next = TListNode_newCopy(nd);
    cur = cur->next;
    object->size ++;
  }
  object->tail = cur;
  return object;
}

TList * TList_appUnique(TList *object, int tag, void * data)
{
  if(!(object->head)){  /* empty list */
    object->head = TListNode_newItem(tag,data);
    object->tail = object->head;
    object->size = 1;
  } else {  /* need to make sure that only unique tags are added */
    TListNode *ndCur = object->head;
    if(ndCur->tag == tag)
      return object;
    while(ndCur->next){
      ndCur = ndCur->next;
      if (ndCur->tag == tag)
	return object;
    }
    object->tail->next = TListNode_newItem(tag,data);
    object->tail= object->tail->next;
    object->size ++;
  }
  return object;
}

TList * TList_append(TList *object, int tag, void *data)
{
  if(!(object->head)){  /* empty list */
    object->head = TListNode_newItem(tag,data);
    object->tail = object->head;
    object->size = 1;
  } else {  
    object->tail->next = TListNode_newItem(tag,data);
    object->tail= object->tail->next;
    object->size ++;
  }
  return object;
}

TList * TList_appTListUnique(TList *object, TList *list)
{
  int fInList;
  
  if(!(object->head)){  /* empty list */
    object=TList_copy(object,list);
  } else {  /* need to make sure that only unique tags are added */
    TListNode *destNode,*cpNode;
    cpNode=list->head;
    while(cpNode){
      destNode=object->head;    /* start at beginning of object list */
      fInList = 0;
      if(destNode->tag == cpNode->tag){
	fInList=1;
      }
      if(!fInList){
	while(destNode->next){
	  destNode = destNode->next;
	  if (destNode->tag == cpNode->tag){
	    fInList=1;
	    break;
	  }
	}
      }
      if(!fInList){                     /* if not in list then add it to end */
	object->tail->next = TListNode_newCopy(cpNode);
	object->tail=object->tail->next;
        object->tail->tag=cpNode->tag;           /* object shouldn't be needed*/
	object->size ++;
      }
      cpNode = cpNode->next;
    }
  }
  return object;
}

TListNode * TList_itemInList(TList *object, int tag)
{
  if(!(object->head)){  /* empty list  */
    return(nullptr);
  } else {  
    TListNode *ndCur = object->head;
    if(ndCur->tag == tag)
      return(ndCur);
    while(ndCur->next){
      ndCur = ndCur->next;
      if (ndCur->tag == tag)
        return(ndCur);
    }
    return(nullptr);
  }
}


TListNode * TList_nodeInList(TList *object, TListNode *node)
{
  if(!(object->head)){  /* empty list  */
    return(nullptr);
  } else {
    TListNode *ndCur = object->head;
    if(ndCur->tag ==node->tag)
      return(ndCur);
    while(ndCur->next){
      ndCur = ndCur->next;
      if (ndCur->tag == node->tag)
        return(ndCur);
    }
    return(nullptr);
  }
}

int TList_size(TList *object) 
    { return object->size; }

int  TList_item(TList *object, int n)
{
TListNode *ndCur = object->head;
int i=0;
  while(i != n){
    i++;
    ndCur = ndCur->next;
  }
  return ndCur->tag;
}

void TList_remove(TList *object, TListNode *rem, TListNode *before)
/* removes the node rem, which is after 'before', if 'before' == NULL
 * then 'rem' is at the head of the list
 */
{
  TListNode *next;

  if(before == nullptr){
    next = object->head->next;
    TListNode_delete(object->head);
    object->head = next;
  } else{
    before->next = rem->next;
    if(object->tail == rem)
      object->tail = before;
    TListNode_delete(rem);
  }
}

void TList_remItem(TList *object, int tag)
/* removes an item from the list
 */
{
  TListNode *last,*ndCur;

  if(!(object->head)){  /* empty list  */
    return;
  } else {  
    last = nullptr;
    ndCur = object->head;
    if(ndCur->tag == tag){
      TList_remove(object,ndCur,last);
      return;
    }
    while(ndCur->next){
      last = ndCur;
      ndCur = ndCur->next;
      if (ndCur->tag == tag){
	TList_remove(object,ndCur,last);
	return;
      }
    }
  }
}
