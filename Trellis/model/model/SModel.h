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
#ifndef H_SModel
#define H_SModel

#include "SModelMember.h"
#include "SSList.h"
#include "SString.h"

namespace AOMD{
class SBoundingBox3d;
}
class SModelMember;
class SModel;

typedef SModel *(*ModelCreatorFunction)(const SString &, int loadOnly, const SSList<SModel*> * const);

/** Base class for all models. At this level a model is just a 
  container for abstract model members. */
class SModel : public SModelMember {
public:
  static SModel * createModel(const SString &name, const SString &type, 
			      const SSList<SModel*> * const = nullptr);
  static void registerModeler(const SString &name,ModelCreatorFunction);

  ~SModel() override;

  /** Get the name of the modeler being used. */
  virtual SString modeler() const =0;
  TopoType::Value type() const override;
  virtual int tag() const;

  SString name() const override;

  AOMD::SBoundingBox3d bounds() const override =0;

  /** Get the model member of the given type with the given tag. */
  virtual SModelMember * getMember(int type, int tag) const = 0;
  /** Get the model member with the given name. If the names are not
   unique, this will return the first one found with this name. */
  virtual SModelMember * getMember(const SString &name) const = 0;
protected:
  SModel(const SString &name);
  SString Name;

  // this is really cheesy, rewrite with stl maps when available
  static int Num;
  static SString Names[10];
  static ModelCreatorFunction IC[10];


};


#endif


