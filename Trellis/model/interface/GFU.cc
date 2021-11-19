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
#include "GFaceUse.h"
#include "GShell.h"
#include "GLoopUse.h"
#include "toPList.h"

extern "C"
pPList GFU_loopUses(GFaceUse *f)
{
  pPList lus = PList_new();
  GLoopUse *lu;
  SSListCIter<GLoopUse*> luIter = f->firstLoopUse();
  while(luIter(lu))
    PList_append(lus,lu);
  return lus;
}

extern "C"
SSListCIter<GLoopUse*> * GFU_loopIter(GFaceUse *lu)
{ return new SSListCIter<GLoopUse*>(lu->firstLoopUse()); }

extern "C"
int GFUIter_next(SSListCIter<GLoopUse*> *luIter, GLoopUse **nlu)
{
  GLoopUse *lu;
  int stat = (*luIter)(lu);
  if(stat){
    *nlu = lu;
    return 1;
  } else 
    delete luIter;
    
  return 0;
}

extern "C"
GRegion * GFU_region(GFaceUse *fu)
{ return fu->shell()->region(); }

extern "C"
GFace * GFU_face(GFaceUse *fu)
{ return fu->face(); }

