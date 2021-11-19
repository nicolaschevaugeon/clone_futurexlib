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
#ifndef _H_TList
#define _H_TList

typedef struct TListNode {
  void * data;
  int tag;
  struct TListNode * next; 
} TListNode;


typedef struct TList{
  TListNode * head;
  TListNode * tail;
  int size;
} TList;

TList * TList_new();

TList * TList_newCopy(TList *);
void TList_delete(TList *);
TList * TList_copy(TList *,TList *source);

TList * TList_appTListUnique(TList *, TList *source);

TList * TList_appUnique(TList *, int tag, void *data);
TList * TList_append(TList *, int tag, void *data);

TListNode * TList_itemInList(TList *, int tag);
TListNode * TList_nodeInList(TList *, TListNode *); 

int TList_size(TList *); 
int TList_item(TList *, int n);

void TList_remove(TList *, TListNode *rem, TListNode *before);

void TList_remItem(TList *, int tag);
#endif
