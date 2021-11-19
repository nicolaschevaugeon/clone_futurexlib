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
#include "GEntity.h"
#include "modeler.h"

void GEN_attachDataP(pGEntity ent, char *tag, void *data) {

  ((AOMD::GEntity *)ent)->setDataP(tag, data);
}

void *GEN_dataP(pGEntity ent, char *tag) {

  return (void *)((AOMD::GEntity *)ent)->dataP(tag);
}

int GEN_modifyDataP(pGEntity ent, char *tag, void *data) {

  return ((AOMD::GEntity *)ent)->modifyDataP(tag, data);
}

void GEN_attachDataI(pGEntity ent, char *tag, int data) {

  ((AOMD::GEntity *)ent)->setDataI(tag, data);
}

int GEN_dataI(pGEntity ent, char *tag) {

  return ((AOMD::GEntity *)ent)->dataI(tag);
}

int GEN_modifyDataI(pGEntity ent, char *tag, int data) {

  return ((AOMD::GEntity *)ent)->modifyDataI(tag, data);
}

void GEN_removeData(pGEntity ent, char *tag) {

  ((AOMD::GEntity *)ent)->removeData(tag);
}
