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

/*
  $Id: ModelDataManager.cc,v 1.9 2005/02/02 00:08:55 acbauer Exp $

File: ModelDataManager.cc
Description: declaration of attached data management class for Model entities
Written by: Jim Teresco
Created: Mon Jun 15 11:50:10 EDT 1998

Modification History:

*/

#include "ModelDataManager.h"
#include "ModelDataId.h"


// unique instance of the ModelDataManager
ModelDataManager modelDataManager;

ModelDataManager::ModelDataManager() : AttachDataManager() {}

ModelDataManager::~ModelDataManager() {

  //AttachDataManager::~AttachDataManager();
}

ModelDataId *ModelDataManager::newModelDataId(const SString &idstr) {

  ModelDataId *new_id=new ModelDataId(idstr);
  registerAttachDataId(new_id);
  return new_id;
}

AttachDataId *ModelDataManager::newAttachDataId(const SString &idstr) {

  return (AttachDataId *)newModelDataId(idstr);
}

void ModelDataManager::deleteModelDataId(ModelDataId *the_id) {

  unregisterId(the_id);
  delete the_id;
}


