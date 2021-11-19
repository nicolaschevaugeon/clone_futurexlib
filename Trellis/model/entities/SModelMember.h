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

#ifndef H_SModelMember
#define H_SModelMember

#include "SBoundingBox3d.h"
#include "AttachableData.h"
#include "ModelDataManager.h"
#include "SSList.h"
#include "GInfo.h"

class GeoRep;
class GRegion;
class GShell;
class GFace;
class GFaceUse;
class GLoopUse;
class GEdge;
class GEdgeUsePair;
class GVertex;
class GVertexUse;
template <class T> class SSList;
class ModelEntityID;
class SString;
class SModel;

/** A member of a model. This may include the model itself or 
one of the entities in the model. An SModelMember is always 
associated with a model (possibly itself)  */
class SModelMember : public AttachableData {
public:
  SModelMember(SModel *model, int tag);
  ~SModelMember() override;

  /** Get the type of this model member. */
  virtual TopoType::Value type() const = 0;

  /** Get the unique, persistent tag associated with this model member. */
  int tag() const;
  int id() const;

  /** Return a string representing the name of the member. The name
    may not be unique */
  virtual SString name() const=0;

  /** Return a bounding box in 3d for the member */
  virtual AOMD::SBoundingBox3d bounds() const =0;

  /** Return true if this contains the given member. */
  virtual int contains(SModelMember *c) const;

  SModel *getModel() const;
  // this is called getModel rather than model so that the derived
  // classes can define model() to return the right type since
  // Sun CC doesn't yet support covarient return types.

  // should only be called by SModelMember operators
  void setID(int id);  
  /* Function to change GTAG such that it corresponds to the same GTAG
     from the Simmetrix model interface. Only used by aomd for sms version 4 file. */
  void setTag(int tag);  
                          

protected:
  
  SModel *Model;
  int ID;
  int d_tag;

  void notImplemented() const;

  AttachDataManager *dataManager() const override { return &modelDataManager; }
};

inline SModel * SModelMember::getModel() const
{ return Model; }

#endif


