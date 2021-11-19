/******************************************************************************** 

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

*********************************************************************************/
#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <map>
#include "AOMD_OwnerManager.h"
#include "mMigrateUtil.h"
#include "ParUtil.h"
#include "mEntity.h"
#include "mMesh.h"


using std::vector;
using std::multimap;
using std::cout;
using std::endl;

namespace AOMD {

AOMD_OwnerManager::AOMD_OwnerManager()
{
  allSharedInfos = new std::multimap <mEntity*,AOMD_SharedInfo>;
  theTopology = new std::multimap <int,int>;
}

AOMD_OwnerManager::~AOMD_OwnerManager()
{
  delete allSharedInfos;
  delete theTopology;
}

void AOMD_OwnerManager::clearRemoteCopies()
{
  allSharedInfos->clear();
}

void AOMD_OwnerManager::addRemoteCopy(mEntity* here, int from, 
                                      mEntity *there , bool per)
{
  allSharedInfos->insert(std::pair<mEntity* const,AOMD_SharedInfo>
                        (here,AOMD_SharedInfo(from,there,per)));
}

void AOMD_OwnerManager::deleteRemoteCopy(mEntity* ent)
{
  iter oit, oitend;
  oit=begin(ent);
  oitend=end(ent);
  allSharedInfos->erase(ent);
}

AOMD_SharedInfo getsecond(const std::pair<mEntity *, AOMD_SharedInfo> & in) {return in.second;};

std::list<AOMD_SharedInfo> AOMD_OwnerManager::getSharedInfo(mEntity * e){
  std::list<AOMD_SharedInfo> listsharedinfo;
  std::transform( begin(e), end(e), std::back_inserter(listsharedinfo), ptr_fun(getsecond) );
  return listsharedinfo;
}

void AOMD_OwnerManager::addNeighbor(int here, int otherside)
{
  theTopology->insert(std::pair<int const,int>(here,otherside));
}

int AOMD_OwnerManager::count(mEntity*e)   
{
  return allSharedInfos->count(e);
}

int AOMD_OwnerManager::count(int i)   
{
  return theTopology->count(i);
}

AOMD_OwnerManager::iter AOMD_OwnerManager:: begin() 
{
  return allSharedInfos->begin();
}

AOMD_OwnerManager::iter AOMD_OwnerManager:: end() 
{
  return allSharedInfos->end();
}

AOMD_OwnerManager::iter_top AOMD_OwnerManager:: begin_top() 
{
  return theTopology->begin();
}

AOMD_OwnerManager::iter_top AOMD_OwnerManager:: end_top() 
{
  return theTopology->end();
}

AOMD_OwnerManager::iter AOMD_OwnerManager:: begin(mEntity*e) 
{
  return allSharedInfos->lower_bound(e);
}

AOMD_OwnerManager::iter AOMD_OwnerManager:: end(mEntity*e) 
{
  return allSharedInfos->upper_bound(e);
}

AOMD_OwnerManager::iter_top AOMD_OwnerManager:: begin_top(int thisside) 
{
  return theTopology->lower_bound(thisside);
}

AOMD_OwnerManager::iter_top AOMD_OwnerManager:: end_top(int thisside) 
{
  return theTopology->upper_bound(thisside);
}

bool AOMD_OwnerManager :: isPeriodic(mEntity *e)
{
  for(iter it = begin(e);it != end(e); ++it)
    if((*it).second.isPeriodic())return true;
  return false;
}

void AOMD_OwnerManager::print()
{
  std::vector<int> ids;
  getPartitionIds (ids);
 
  for(int i=0;i<ids.size();++i)
    {
      printf("Processor %d boundary %d comm (",ParUtil::Instance()->rank(),ids[i]);
      for(iter_top itx = begin_top(ids[i]); itx != end_top(ids[i]); ++itx)
	{	      
	  int to = (*itx).second;
	  printf("%d ",to);
	}
      printf(")\n");
    }
//  return;
  for(iter it = begin() ; it != end() ; ++it)
    {
      AOMD_SharedInfo si = (*it).second;
      printf("connect with P%d",si.pid()); 
      mEntity *e = (*it).first;
      e->print();
    }

}


void AOMD_OwnerManager::getPartitionIds ( std::vector<int> &ids )
{
  int id_old = 0;
  for(iter_top itx = begin_top(); itx != end_top(); ++itx)
    {	      
      int id = (*itx).first;
      if(id != id_old)ids.push_back(id);
      id_old = id;
    }
}

int AOMD_OwnerManager::find (int nbProcs , int *procs)
{

  std::vector<int> ids;
  std::sort(procs,procs+nbProcs);
  getPartitionIds(ids);
	
  // MUST ADD THIS PROCESSOR TO THE LIST

  for(int i=0;i<ids.size();++i)
    {
      if(nbProcs == count(ids[i])) // if size corresponds 
	{
	  bool found = true;
	  for(iter_top itx = begin_top(ids[i]); itx != end_top(ids[i]); ++itx)
	    {	      
	      int to = (*itx).second;
	      if(!std::binary_search(procs,procs+nbProcs,to))found = false;	      
	    }
	  if(found) return ids[i];
	}
    }
  return 0;
}

int AOMD_OwnerManager::create (int nbProcs , int *procs)
{
  int maxid = 0;
  for(iter_top itx = begin_top(); itx != end_top(); ++itx)
    {	      
      int id = (*itx).first;
      if(id > maxid)maxid = id;
    }
  for(int i=0;i<nbProcs;i++)
    {
      addNeighbor(maxid+1,procs[i]);
    }
  return maxid+1;
}


void AOMD_OwnerManager::setPoorPid_asOwner(mMesh* mesh)
{

}
  
bool AOMD_OwnerManager::amITheOwnerOf (mEntity * e )
{
  if (e->getOwner()==ParUtil::Instance()->rank())
    return true;
  else
    return false;
}

} // end of namespace
