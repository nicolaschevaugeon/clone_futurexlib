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
#include "ModelEntityID.h"
#include "SModelMember.h"
#include "SModel.h"


#include <iostream>
using namespace std;


#include <cassert>

#ifdef OTL
#include "otl.h"
#include "DBconnection.h"
#endif

ModelEntityID::ModelEntityID(SModelMember *ent)
{
  Type = ent->type();
  Tag = ent->tag();
}

#ifdef OTL
ModelEntityID::ModelEntityID(int ModelAssocID)
{
  DBconnection *dbc = DBconnection::Instance();
  try {
    otl_stream i(50,"select ModelEntityType, ModelEntityTag from ModelEntityAssociation where ModelAssocId=:f<int>",*(dbc->getHandler()));
    i << ModelAssocID;
    int counter = 0;
    while (!i.eof())
    {
      counter++;
      i >> Type >> Tag;
    }
    assert(counter==1);
  }
  catch(otl_exception &p)
  {
    std::cerr << p.msg<< std::endl;
    std::cerr<<p.stm_text<< std::endl;
    std::cerr << p.var_info << std::endl;
  }
}
#endif

ModelEntityID::ModelEntityID(istream &in)
{
  in >> Type >> Tag;
  if(Type < 0)
    cerr << "ModelEntityID::ModelEntityID(istream &in) - error, invalid type\n";
}

int ModelEntityID::operator==(const ModelEntityID &e) const
{
  return (Type == e.Type && Tag == e.Tag);
}

int ModelEntityID::operator!=(const ModelEntityID &e) const
{
  return !operator==(e);
}

SModelMember *ModelEntityID::entity(SModel *model) const
{
  return model->getMember(Type,Tag);
}

void ModelEntityID::write(ostream &os) const
{
  os << Type << " " << Tag;
}

#ifdef OTL
void ModelEntityID::writeToDB(int ModelAssocID) const
{
  DBconnection *dbc = DBconnection::Instance();
  try {
    otl_stream o(50,"insert into ModelEntityAssociation values( :f1<int>, :f2<int>, :f3<int>)",*(dbc->getHandler()));
    o << ModelAssocID << Type << Tag;
  }
  catch(otl_exception &p)
  {
    std::cerr << p.msg<< std::endl;
    std::cerr<<p.stm_text<< std::endl;
    std::cerr << p.var_info << std::endl;
  }
}
#endif


