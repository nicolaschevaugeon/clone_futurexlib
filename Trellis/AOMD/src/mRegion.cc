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
#include <iostream>
#include "mRegion.h"
#include "mEdge.h"
#include "mException.h"

namespace AOMD {

mEdge * mRegion::commonEdge (mRegion *r1, mRegion *r2)
{
  std::set<mEdge*,EntityLessThanKey> thisEdges;
  std::set<mEdge*,EntityLessThanKey> r1Edges;
  std::set<mEdge*,EntityLessThanKey> r2Edges;
 
  for(int i=0;i<size(1);i++)thisEdges.insert((mEdge*)get(1,i));
{  for(int i=0;i<r1->size(1);i++)r1Edges.insert((mEdge*)r1->get(1,i));}
{  for(int i=0;i<r2->size(1);i++)r2Edges.insert((mEdge*)r2->get(1,i));}

  for( std::set<mEdge*,EntityLessThanKey>::const_iterator iter = thisEdges.begin();
       iter != thisEdges.end();
       ++iter)
    {
      if(r1Edges.find(*iter) != r1Edges.end() &&
	 r2Edges.find(*iter) != r2Edges.end())return *iter;
    }
  return nullptr;


  throw new mException (__LINE__, __FILE__,"mRegion::commonEdge not none yet ...");
}

mFace * mRegion::commonFace (mRegion *r)
{
  for(int i=0;i<size(2);i++)
    {
      mEntity *f = get(2,i);
      for(int j=0;j<r->size(2);j++)
	{
	  if(f == r->get(2,j))return (mFace*)f;
	}      
    }
  return nullptr;
}

} // end of namespace
