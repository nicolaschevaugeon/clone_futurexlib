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
  $Id: ModelDataManager.h,v 1.7 2005/02/02 00:08:55 acbauer Exp $

File: ModelDataManager.h
Description: Header file for attached data management class for mesh entities
Written by: Jim Teresco
Created: Mon Jun 15 11:50:10 EDT 1998

Modification History:

*/

#ifndef H_ModelDataManager
#define H_ModelDataManager

//#include "modelfwd.h"
//#include "utilfwd.h"
#include "AttachDataManager.h"

#include "ModelDataId.h"


class ModelDataManager: public AttachDataManager {
public:
  ModelDataManager();
  ~ModelDataManager() override;

  // virtual member to return an DataId as the base AttachDataId
  AttachDataId *newAttachDataId(const SString &) override;
  virtual ModelDataId *newModelDataId(const SString &);
  void deleteModelDataId(ModelDataId *);
};

// reference to the unique instance of the ModelDataManager
extern ModelDataManager modelDataManager;

#endif
