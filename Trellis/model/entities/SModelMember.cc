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
#include "SModelMember.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
//#include "GShell.h"
#include "GRegion.h"
#include "SSList.h"
#include "SGModel.h"
#include "MessageOut.h"

SModelMember::SModelMember(SModel *model, int t)
  : Model(model), d_tag(t)
{ }

SModelMember::~SModelMember() 
= default;

int SModelMember::tag() const
{ return d_tag; }

int SModelMember::id() const
{ return ID; }

void SModelMember::setID(int id)
{ ID = id; }

void SModelMember::setTag(int tag)
{ d_tag = tag; }

void SModelMember::notImplemented() const
{ Warning("Operation not implemented"); }

int SModelMember::contains(SModelMember *c) const
{ return this==c; }






