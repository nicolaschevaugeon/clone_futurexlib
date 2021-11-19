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
#ifndef H_ListMan
#define H_ListMan

/* $Id: ListMan.h,v 1.3 2005/02/02 00:08:57 acbauer Exp $ */
#include "PList.h"

typedef struct ListMan{
  int Num;
  int MaxNum;
  PList **lists;
} ListMan;

#ifdef __cplusplus
extern "C" {
#endif

void ListMan_init(int max);

void ListMan_quit(void);

PList * ListMan_get(void);
void ListMan_put(PList *list);

#ifdef __cplusplus
}
#endif

#endif
