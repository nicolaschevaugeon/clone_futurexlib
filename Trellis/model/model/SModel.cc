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
#include "SModel.h"
#include <cstdlib>


using std::cerr;

int SModel::Num = 0;
SString SModel::Names[10];
ModelCreatorFunction SModel::IC[10];

SModel * SModel::createModel(const SString &name, const SString &type,
			     const SSList<SModel*> * const deps)
{
  for(int i = 0; i < Num; i++){
    if(type == Names[i])
      return IC[i](name,0,deps);
  }
  cerr << "Can't create model of type " << type << "\n";
  return nullptr;
}

void SModel::registerModeler(const SString &name,ModelCreatorFunction f)
{
  if(f){
    Names[Num] = name;
    IC[Num] = f;
    Num++;
    if(Num > 10){
      cerr << "SModel dictionary overflow\n";
      exit(0);
    }
  }
}


SModel::SModel(const SString &name)
: SModelMember(nullptr,0), Name(name)
{
}

SModel::~SModel()
= default;

SString SModel::name() const
{ return Name; }

TopoType::Value SModel::type() const
{ return TopoType::Model; }

int SModel::tag() const
{ return 0; }
