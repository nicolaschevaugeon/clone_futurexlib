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
 *  AOMD/include/pmEntity.h
 *  Created by Seegyoung Seol, on Mon Dec 08 2003, 10:26:37 EDT
 *
 *  File Content: partition model entity class
 *
 *************************************************************************** </i>*/

#ifndef _PM_ENTITY_H_
#define _PM_ENTITY_H_

#include <set>
#include <vector>

namespace AOMD {
  
  class pmEntity;
  
  class pmEntityLessThanKey
  {
  public :
    bool operator() (pmEntity *cb1, pmEntity* cb2) const;
  };
  
  class pmEntity
  {
private:
  int iD;
  int level;	// dimension
  int owner;		// owner partition
  std::set<int> BPs;	// bounding partitions

public:
  typedef std::set<int>::iterator BPIter; 

// constructor + destructor  
  pmEntity(int i, int l): iD(i), level(l) {}
  pmEntity(int i, std::set<int>& bps, int l);
//  ~pmEntity();

// get data members
  int getId() const {return iD;}
  int getLevel() const {return level;}
  int getOwner() { return owner; }

// set data members
  void setOwner(int o) { owner=o; }
  void setOwner(std::vector<int>&);

// equality check based on BPS
  bool isEqual(std::set<int>& bpSet);
  void getBPs(std::vector<int>&);
  void addBPs(int p) {  BPs.insert(p); }    
  int getNumBPs() { return BPs.size(); }

// BPs iterator
  
  BPIter bpBegin() { return BPs.begin();  }
  BPIter bpEnd()   { return BPs.end();  }

// accessory
  void print();
  
  /*  protected:
    int iD;
    int level;
    int owner;
  public:
    pmEntity(int i, int l): iD(i), level(l) {}
    //~mCommonBdry();
    int getId() const {return iD;}
    int getLevel() const {return level;}
    void setOwner(int o) { owner=o; }
    int getOwner() { return owner; }
    */
  };  
} // end of namespace

#endif
