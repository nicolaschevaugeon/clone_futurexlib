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
/* imported from SCOREC Mesh Database, July 1997 */
/* $Id: PList.cc,v 1.3 2005/02/02 00:08:58 acbauer Exp $ */

#include "PList.h"

#include "ListMan.h"
#include "MD_memory.h"
#include <cstdio>
#include <cstdlib>


PList * PList_allocate(void)
{ /* OK */
  PList *cthis;
  cthis = (PList *)MD_malloc(sizeof(PList));
  cthis->ent = (void**)MD_malloc(4*sizeof(void *));
  cthis->size = 0;
  cthis->allocSize = 4;
  return cthis;
} 

PList * PList_new(void)
{ /* OK */
  return ListMan_get();
}

void PList_growToSize(PList *l, int size)
{
  int i;
  if(l->allocSize >= size)
    return;
  else{
    void **elist = l->ent;
    l->ent =  (void**)MD_malloc(size*sizeof(void *));
    for(i=0; i < l->size; i++)
      l->ent[i] = elist[i];
    MD_free(l->allocSize*sizeof(void*),elist);
    l->allocSize = size;
  }
}
  
PList * PList_newCopy(PList *il)
{ /* OK */
  int i;
  PList *cthis;
  void **e;
  void **ec;
  cthis = ListMan_get();

  if(il->size){
    PList_growToSize(cthis,il->size);
    e = il->ent;
    ec = cthis->ent;
    for(i=0; i < il->size; i++)
      *ec++ = *e++;
  }
  cthis->size = il->size;
  return cthis;
}

void PList_deallocate(PList *cthis)
{ 
  MD_free(cthis->allocSize*sizeof(void*),cthis->ent);
  MD_free(sizeof(PList),cthis);
}

void PList_delete(PList *cthis)
{ /* OK */
  PList_clear(cthis);
  ListMan_put(cthis);
}

void PList_clear(PList *cthis)
{ /* OK */
  cthis->size=0;
}

PList * PList_copy(PList *cthis,PList *il)
{ /* OK */
  int i;
  void **e, **ec;

  if(il->size){
    PList_growToSize(cthis,il->size);
    e = il->ent;
    ec = cthis->ent;
    for(i=0; i < il->size; i++)
      *ec++ = *e++;
    cthis->size = il->size;
  }
  return cthis;

}

PList * PList_fill(PList *cthis,void **array, int n)
/* fills a Plist from an array of size n */
{ 
  int i = 0;
  PList_growToSize(cthis,n);
  for(i=0; i < n; i++)
    PList_append(cthis,array[i]);
  return cthis;
}

PList * PList_appUnique(PList *cthis, void *tag)
{ /* OK */
  if( ! itemInList(cthis,tag))
    cthis = PList_append(cthis,tag);
  return cthis;
}

PList * PList_append(PList *cthis, void *tag)
{ /* OK */
  
  if(cthis->size+1 > cthis->allocSize)
    PList_growToSize(cthis,2*(cthis->size)+1); /* +1 in case size==0 */
  
  cthis->ent[cthis->size] = tag;  
  cthis->size++;
  return cthis;
}

PList * PList_appPList(PList *cthis, PList *tagL)
{ /* OK */
  int i;

  for(i=0; i < PList_size(tagL); i++)
    PList_append(cthis,PList_item(tagL,i));
  return cthis;
}

PList * PList_appPListUnique(PList *cthis, PList *tagL)
{ /* OK */
  int i;

  for(i=0; i < PList_size(tagL); i++)
    PList_appUnique(cthis,PList_item(tagL,i));
  return cthis;
}

int itemInList(PList *cthis, void *tag)
{ /* OK */
  int i;
  void **e;
  e = cthis->ent;
  for(i=0; i < cthis->size; i++){
    if( *e++ == tag)
      return 1;
  }
  return 0;
}

int PList_inList(PList *cthis, void *tag)
{ /* OK */
  int i;
  void **e;
  e = cthis->ent;
  for(i=0; i < cthis->size; i++){
    if( *e++ == tag)
      return 1;
  }
  return 0;
}

int PList_replace(PList *cthis, void *tag, void *newtag)
{ /* OK */
  int i;
  void **e;
  e = cthis->ent;
  for(i=0; i < cthis->size; i++){
    if( *e == tag){
      *e = newtag;
      return 1;
    }
    e++;
  }
  return 0;
}

int PList_replaceNth(PList *cthis, void *newtag, int n)
{
  cthis->ent[n] = newtag;
  return 1;
}

int PList_size(PList *cthis)  /* OK */
    { return cthis->size; }

void *PList_next(PList *cthis , void **temp)
{
  /*void *item;*/

/*
  if(*temp == 0){
    *temp = malloc(sizeof(int));
    *(int*)(*temp) = 0;
  }
*/
  if(*(int*)(temp) >= cthis->size){
    /*free(*temp);*/
    return nullptr;
  }
  return cthis->ent[(*(int*)(temp))++];
}

void *PList_item(PList *cthis, int n)
{ /* OK */
  return cthis->ent[n];
}

void PList_remItem(PList *cthis, void *tag)
{ /* OK */
  int i;
  int move = 1;
  for(i=0; i < cthis->size; i++)
    if(cthis->ent[i] == tag){
      cthis->size--;
      break;
    }
  for( ; i < cthis->size; i++){
    if(cthis->ent[i+1] == tag){
      move++;
      cthis->size--;
      continue;
    }
    cthis->ent[i] = cthis->ent[i+move];
  }
}

void PList_remFirst(PList *cthis, void *tag)
{ /* OK */
  int i;
  for(i=0; i < cthis->size; i++)
    if(cthis->ent[i] == tag){
      cthis->size--;
      break;
    }
  for( ; i < cthis->size; i++)
    cthis->ent[i] = cthis->ent[i+1];
}


void PList_forEach(PList *cthis, void func(void *) )
{
  int i;
  for(i=0; i < cthis->size; i++)
    func(cthis->ent[i]);
}

#if 0
#ifdef FORT_LC_US
void plist_new_(PList **list)
#endif
#ifdef FORT_LC
void plist_new(PList **list)
#endif
{
  *list = PList_new();
}

#ifdef FORT_LC_US
void plist_delete_(PList **list)
#endif
#ifdef FORT_LC
void plist_delete(PList **list)
#endif
{
  PList_delete(*list);
}

#ifdef FORT_LC_US
void plist_size_(PList **list, int *s)
#endif
#ifdef FORT_LC
void plist_size(PList **list, int *s)
#endif
{
  *s = PList_size(*list);
}

#ifdef FORT_LC_US
void plist_item_(PList **cthis, int *n, Entity **e)
#endif
#ifdef FORT_LC
void plist_item(PList **cthis, int *n, Entity **e)
#endif
{
  *e = PList_item(*cthis,*n);
}
#endif

