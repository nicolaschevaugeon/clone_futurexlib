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
/* $Id: ListMan.cc,v 1.2 2005/02/02 00:08:57 acbauer Exp $ */
#include "ListMan.h"
#include "PList.h"
#include <cstdio>
#include <cstdlib>

ListMan gListMan;  /* the global list manager object, do we need one for 
		      each type ?? */

void MD_debugListCacheOff()
{
  gListMan.MaxNum = 0;
}

void ListMan_init(int max)
{
  gListMan.MaxNum = max; 
  gListMan.Num = 0;
  gListMan.lists = (PList **)malloc(max * sizeof(PList*));
}

void ListMan_quit(void)
{
  int i;

  for(i=0; i < gListMan.Num; i++){
    if(gListMan.lists[i])
      PList_deallocate(gListMan.lists[i]);
  }
  free(gListMan.lists);
  gListMan.lists = nullptr;
  gListMan.MaxNum = 0; 
  gListMan.Num = 0;
  
}

PList * ListMan_get(void)
{
  PList *list;

  if(gListMan.Num){
    list = gListMan.lists[gListMan.Num-1];
    gListMan.lists[gListMan.Num-1] = nullptr;
    gListMan.Num --;
  } else
    list = PList_allocate();

  return list;
}
    

void ListMan_put(PList *list)
{
  if(gListMan.Num < gListMan.MaxNum){
    gListMan.Num++;
    gListMan.lists[gListMan.Num-1] = list;
  } else
    PList_deallocate(list);
}
