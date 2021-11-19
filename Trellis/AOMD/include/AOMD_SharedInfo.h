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
#ifndef _AOMD_SHAREDINFO_H_
#define _AOMD_SHAREDINFO_H_

namespace AOMD {
  
  class mEntity;
  
  struct ltshared
  {
    bool operator()(const mEntity* e1, const mEntity* e2) const
    {
      return e1 < e2;
    }
  };
  
  class AOMD_SharedInfo
  {
  private:
    int PID;
    mEntity *remote;
    bool periodic;
  public:
    inline AOMD_SharedInfo(int p, mEntity* e, bool per)
      :PID(p),remote(e), periodic(per)
    {}
    inline mEntity *getRemotePointer() {return remote;}
    inline mEntity *getRemotePointer()const  {return remote;}

    inline int pid() const {return PID;}
    inline bool isPeriodic() const {return periodic;}
  };
}
#endif
