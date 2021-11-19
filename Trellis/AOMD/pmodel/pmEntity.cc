/**************************************************************************** 

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written 
   and maintained by the Scientific Computation Research Center (SCOREC) at 
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*****************************************************************************/

/*<i> ****************************************************************************
 *
 *  AOMD/pmodel/pmEntity.cc
 *  Created by Seegyoung Seol, on Mon Dec 08 2003, 10:26:37 EDT
 *
 *  File Content: partition model entity class definition
 *
 *************************************************************************** </i>*/

#include "pmEntity.h"
#include "ParUtil.h"
#include <iostream>
#include <set>
#include <algorithm>
#include <vector>
#include <cassert>

using std::set;
using std::cout;
using std::endl;
using std::vector;

namespace AOMD {

bool pmEntityLessThanKey::operator() (pmEntity* cb1, pmEntity* cb2) const
{
  if(cb1->getId() < cb2->getId())return true;
  if(cb1->getId() >= cb2->getId())return false;
//    if(cb1->getLevel() < cb2->getLevel())return true;
//    if(cb1->getLevel() > cb2->getLevel())return false;

  return false;
}

pmEntity::pmEntity(int id, set<int>& bps, int dim)
{
  iD = id;
  level = dim;
  owner=-1;
  for (set<int>::iterator it=bps.begin(); it != bps.end(); ++it)
    BPs.insert(*it);
}

bool pmEntity::isEqual(set<int>& bpSet)
{
  if (BPs.size()!=bpSet.size())
    return false;
  set<int>::iterator it1 = BPs.begin();
  set<int>::iterator it2 = bpSet.begin();
  for (; it1!=BPs.end();++it1, ++it2)
  {
    if ((*it1) != (*it2))
      return false;
  }
  return true;
//  return !(lexicographical_compare(begin(), end(),
// 				  bpSet.begin(), bpSet.end()));
}

void pmEntity::setOwner(std::vector<int>& allSortedPids)
{
  vector<int>::iterator vit;
  set<int>::iterator where;
  for (vit=allSortedPids.begin(); vit!=allSortedPids.end();++vit)
  {
    where = BPs.find(*vit);
    if (where==BPs.end())
      continue;
    owner = *vit;
    return;
  }
  assert(false);  // if reach this line, error!
}

void pmEntity::print()
{
  cout<<"("<<ParUtil::Instance()->rank()<<") PE iD="<<iD<<", BPs=(";
  for (BPIter it=bpBegin(); it!=bpEnd(); ++it)
    cout<<*it<<",";
  cout<<"), dim="<<level<<", owner="<<owner<<endl;
}
  
void pmEntity::getBPs(vector<int>& bps)
{
  std::copy(bpBegin(),bpEnd(),back_inserter(bps));
}

} // end of namespace
