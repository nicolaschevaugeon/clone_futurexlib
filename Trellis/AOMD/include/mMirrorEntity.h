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
#ifndef _MMIRRORVERTEX_H_
#define _MMIRRORVERTEX_H_

/*
  Entity which is like a mirror for other ones
*/
#include <vector>
#include "mVertex.h"

namespace AOMD {

class mVertex;

class mMirrorVertex 
//nico : public mVertex
{
protected:
  std::vector<mVertex*> copies;
public:
  // for mirror Vertices
  //nico mMirrorVertex(int iD);
  virtual ~mMirrorVertex()= default;
  void addCopy (mVertex *);
  inline int nbCopies () const {return copies.size();} 
  inline mVertex* getCopy  (int i) {return copies[i];} 
	// ************* TEMPORARY : COMPATIBILITY WITH GENMAX *********************
	inline void Clear() { copies.clear(); } 
	// ************** CENAERO (OM-EW) ******************************************
  //nico virtual mEntity::mType getType() const{return MIRROR;}
  virtual void print() const;
};

} // end of namespace

#endif 
