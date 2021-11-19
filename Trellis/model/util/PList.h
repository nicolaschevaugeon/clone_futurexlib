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
/* $Id: PList.h,v 1.6 2005/02/02 00:08:58 acbauer Exp $ */
#ifndef _H_PList
#define _H_PList

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PList{
  void **ent;
  int size;
  int allocSize;
} PList;

typedef struct PList * pPList;

PList * PList_new(void);

PList * PList_newCopy(PList *);
void PList_delete(PList *);
PList * PList_copy(PList *,PList *source);
PList * PList_fill(PList *this_one,void **array, int n);

PList * PList_appPListUnique(PList *, PList *source);
PList * PList_appPList(PList *, PList *source);

PList * PList_appUnique(PList *, void *item);
PList * PList_append(PList *, void *item);

int PList_inList(PList *, void *tag);
int itemInList(PList *, void *);
int PList_replace(PList *, void *tag, void *newtag);
int PList_replaceNth(PList *this_one, void *newtag, int n);
void PList_remCurrent(PList *, void **temp);

int PList_size(PList *); 
void * PList_item(PList *, int n);
void * PList_next(PList *, void **temp);

void PList_remItem(PList *, void *item);
void PList_remFirst(PList *, void *tag);

void PList_forEach(PList *, void func(void *) );

void PList_clear(PList *);
PList * PList_allocate(void);
void PList_deallocate(PList *);
void PList_print(PList *this_one);

#ifdef __cplusplus
}
#endif

#endif
