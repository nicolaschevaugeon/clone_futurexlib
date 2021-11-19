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

/*****************************************************************************
 *
 *  include/mCommonBdry.h
 *  Created by Eunyoung Seol, on Mon Aug 25 2003, 10:26:37 EDT
 *
 *  File Content: common boundary class
 ****************************************************************************/

#ifndef _MCOMMON_BDRY_H_
#define _MCOMMON_BDRY_H_

namespace AOMD {
  
  class mCommonBdry;
  
  class mCommonBdryLessThanKey
  {
  public :
    bool operator() (mCommonBdry *cb1, mCommonBdry* cb2) const;
  };
  
  class mCommonBdry
  {
  protected:
    int iD;
    int level;
    int owner;
  public:
    mCommonBdry(int i, int l): iD(i), level(l) {}
    //~mCommonBdry();
    int getId() const {return iD;}
    int getLevel() const {return level;}
    void setOwner(int o) { owner=o; }
    int getOwner() { return owner; }
  };  
} // end of namespace

#endif
