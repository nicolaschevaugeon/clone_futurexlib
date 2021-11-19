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
// mPrism.cpp: implementation of the mPrism class.
#include <iostream>
#include "mPrism.h"
#include "mEntityContainer.h"
#include "mException.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mFace.h"

#include <cstdio>

namespace AOMD {
  
  mPrism::mPrism()
    : mRegion()
  {
  }
  
  mPrism::mPrism(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4,
		 mVertex *v5, mVertex *v6, 
		 GEntity *classification)
    : mRegion()
  {
    theAdjacencies[0] = new mAdjacencyContainer (6);
    theAdjacencies[0]->add(v1);
    theAdjacencies[0]->add(v2);	
    theAdjacencies[0]->add(v3);	
    theAdjacencies[0]->add(v4);	
    theAdjacencies[0]->add(v5);	
    theAdjacencies[0]->add(v6);	
    theClassification = classification;
    iD = v1->getRAND() + v2->getRAND() + v3->getRAND() + v4->getRAND() +
    v5->getRAND() + v6->getRAND();
  }
  
  int mPrism:: getNbTemplates (int what) const
  {
    switch(what)
      {
      case 0: return 6;
      case 1: return 9;
      case 2: return 5;
	
      default : 
	{
	  char text[256];
	  sprintf(text,"trying to ask for adjacency list of level %d for a prism when it's not created",what);
	  throw new mException (__LINE__,__FILE__,text);
	}
      }
  }

  static int Tev[9][2] = {{0,1},{0,2},{1,2},
			  {3,4},{3,5},{4,5},
			  {0,3},{1,4},{2,5}};

  static int Tfv[5][4] =   {{0,1,2,-1},
			    {3,5,4,-1},
			    {0,3,4,1},
			    {1,4,5,2},
			    {0,2,5,3}};

  static int Tfe[5][4] =   {{0,1,2,-1},
			    {3,4,5,-1},
			    {6,3,7,0},
			    {7,5,8,2},
			    {6,1,8,4}};

  static int Tve[6][2] = {{0,1},{0,2},{1,2},{0,3},{3,5},{4,5}};
  static int Tvf[6][3] = {{0,2,4},{0,2,3},{0,3,4},{1,2,4},{1,2,3},{1,3,4}};
  static int Tef[9][2] = {{0,2},{0,4},{0,3},
			   {1,2},{1,4},{1,3},
			   {2,4},{2,3},{3,4}};
  
  mEntity *mPrism::getTemplate(int ith, int what, int with)const
  {
    switch(what)
      {
      case 0:
	if(theAdjacencies[0])return (*theAdjacencies[0])[ith];
	if(theAdjacencies[1])
	  {
	    mEdge *e1 = (mEdge*)get(1,Tve[ith][0]);
	    mEdge *e2 = (mEdge*)get(1,Tve[ith][1]);
	    return e1->commonVertex(e2);
	  }
	else if (theAdjacencies[2])
	  {
	    mFace *f1 = (mFace*)get(2,Tvf[ith][0]);
	    mFace *f2 = (mFace*)get(2,Tvf[ith][1]);
	    mFace *f3 = (mFace*)get(2,Tvf[ith][2]);
	    return f1->commonVertex(f2,f3);
	  }
	break;
      case 1:
	if (theAdjacencies[2])
	  {
	    mFace *f1 = (mFace*)get(2,Tef[ith][0]);
	    mFace *f2 = (mFace*)get(2,Tef[ith][1]);
	    mEdge *e =  f1->commonEdge(f2);
	    if(e) return e;
	  }
	return new mEdge ((mVertex*)theAdjacencies[0]->get(Tev[ith][0]),
			  (mVertex*)theAdjacencies[0]->get(Tev[ith][1]),
			  theClassification);
	break;
      case 2:
	if(Tfv[ith][3] == -1)
	  {
	    if(with == 0)
	      return new mFace ((mVertex*)theAdjacencies[0]->get(Tfv[ith][0]),
				(mVertex*)theAdjacencies[0]->get(Tfv[ith][1]),
				(mVertex*)theAdjacencies[0]->get(Tfv[ith][2]),
				theClassification);
	    else if (with == 1)
	      {
		return new mFace ((mEdge*)theAdjacencies[1]->get(Tfe[ith][0]),
				  (mEdge*)theAdjacencies[1]->get(Tfe[ith][1]),
				  (mEdge*)theAdjacencies[1]->get(Tfe[ith][2]),
				  theClassification);
	      }
	    else throw new mException (__LINE__,__FILE__,"weird adjacency asked");
	  }
	else
	  {
	    if(with == 0)
	      return new mFace ((mVertex*)theAdjacencies[0]->get(Tfv[ith][0]),
				(mVertex*)theAdjacencies[0]->get(Tfv[ith][1]),
				(mVertex*)theAdjacencies[0]->get(Tfv[ith][2]),
				(mVertex*)theAdjacencies[0]->get(Tfv[ith][3]),
				theClassification);
	    else if (with == 1)
	      {
		return new mFace ((mEdge*)theAdjacencies[1]->get(Tfe[ith][0]),
				  (mEdge*)theAdjacencies[1]->get(Tfe[ith][1]),
				  (mEdge*)theAdjacencies[1]->get(Tfe[ith][2]),
				  (mEdge*)theAdjacencies[1]->get(Tfe[ith][3]),
				  theClassification);
	      }
	    else throw new mException (__LINE__,__FILE__,"weird adjacency asked");
	  }
      }
    return nullptr;
  }
  
  mPrism::~mPrism()
  = default;
  
} // end of namespace
