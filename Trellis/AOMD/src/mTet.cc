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
#include "AOMDfwd.h"
#include "AOMD_Defs.h"
#include "mTet.h"
#include "mEntityContainer.h"
#include "mException.h"
#include "mVertex.h"
#include "mVector.h"
#include "mEdge.h"
#include "mFace.h"
#include "AOMD_Internals.h"
#include <cstdio>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

namespace AOMD {

  static double Volume (mTet *tet)
  {
    mVertex *v1 = (mVertex*)tet->get(0,0);
    mVertex *v2 = (mVertex*)tet->get(0,1);
    mVertex *v3 = (mVertex*)tet->get(0,2);
    mVertex *v4 = (mVertex*)tet->get(0,3);

    Trellis_Util::mVector v12 (v2->point(),v1->point());
    Trellis_Util::mVector v13 (v3->point(),v1->point());
    Trellis_Util::mVector v14 (v4->point(),v1->point());
    
    return -(1./6.) * ((v12 % v13) * v14);
  }
  
  mTet::mTet()
    : mRegion()
  {
  }
  
  mTet::mTet(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, GEntity *classification)
    : mRegion()
  {
    theAdjacencies[0] = new mAdjacencyContainer (4);
    theAdjacencies[0]->add(v1);
    theAdjacencies[0]->add(v2);	
    theAdjacencies[0]->add(v3);	
    theAdjacencies[0]->add(v4);	
    theClassification = classification;
    iD = v1->getRAND() + v2->getRAND() + v3->getRAND() + v4->getRAND();
    //uiD = buildUid();
    //iD = v1->getId() + v2->getId() + v3->getId() + v4->getId();
    double volume = Volume(this);
//    if (volume < 1.e-12 && volume>0.0)
//      {
//	v1->print();
//	v2->print();
//	v3->print();
//	v4->print();
       // cout<<"("<<M_Pid()<<") WARNING: sliver region found\n    ";
       // print(); // meaning sliver region
      //cout<<"("<<M_Pid()<<") WARNING: sliver region "<<getUid()<<" created\n    ";
//      }

    if(volume < 0.0)
      {
	delete theAdjacencies[0];
	theAdjacencies[0] = new mAdjacencyContainer (4);
	theAdjacencies[0]->add(v2);
	theAdjacencies[0]->add(v1);	
	theAdjacencies[0]->add(v3);	
	theAdjacencies[0]->add(v4);	
        //cout<<"("<<M_Pid()<<") ERROR: negative region "<<getUid()<<" created\n" ;      
      } 
  }

mTet::mTet(mFace *f1, mFace *f2, mFace *f3, mFace *f4, GEntity *classification)
  : mRegion ()
{

    // check & volume sign => need 4 tet vertex
    mVertex * v1 = f1->commonVertex(f2,f4);
    mVertex * v2 = f1->commonVertex(f2,f3);
    mVertex * v3 = f1->commonVertex(f3,f4);
    mVertex * v4 = f2->commonVertex(f4,f3);
  
    // check all face are well connected 3 by 3
    if(!v1 || !v2 || !v3 || !v4 ) { throw new mException (__LINE__,__FILE__,"3 faces does not have a common vertex"); }

    // vector to compute sign
    Trellis_Util::mVector v12 (v2->point(),v1->point());
    Trellis_Util::mVector v13 (v3->point(),v1->point());
    Trellis_Util::mVector v14 (v4->point(),v1->point());
    
    theAdjacencies[2] = new mAdjacencyContainer (4);
    // depending on sign  adjancy is modifyed from input. Dangerous ! But it was like that with a bugg at time of commit. Consistant 
    // with constructor with vertex.
    if ( ((v12 % v13) * v14)> 0.0 )
    {
        theAdjacencies[2]->add(f2);
        theAdjacencies[2]->add(f1);	
    }
    else
    {
        theAdjacencies[2]->add(f1);
        theAdjacencies[2]->add(f2);	
    }
    theAdjacencies[2]->add(f3);	
    theAdjacencies[2]->add(f4);	
    theClassification = classification;
    iD = v1->getRAND() + v2->getRAND() + v3->getRAND() + v4->getRAND();

}

int mTet:: getNbTemplates (int what) const
{
  switch(what)
    {
    case 0: return 4;
    case 1: return 6;
    case 2: return 4;
    default : throw new mException (__LINE__,__FILE__,"");
    }
}

// --- DIRECT TABLES ---
// gives edges as a function of vertices
// static int Tev[6][2] = {{0,1},{0,2},{0,3},{1,2},
//			{1,3},{2,3}};
static int Tev[6][2] = {{0,1},{1,2},{2,0},{0,3},
			{1,3},{2,3}};

// gives faces as a function of vertices
//  static int Tfv[4][3] =   {{0,1,2},
//  			  {0,1,3},
//  			  {0,2,3},
//  			  {1,2,3}};
static int Tfv[4][3] =   {{0,1,2},
			  {0,1,3},
			  {1,2,3},
			  {0,2,3}};

// gives faces as a function of edges (composition 
// of the 2 first ones)
// static int Tfe[4][3] = {{0,3,1},{0,4,2},{1,5,2},{3,5,4}};
static int Tfe[4][3] = {{0,1,2},{0,4,3},{1,5,4},{2,3,5}};

// --- INVERT TABLES ---
// gives vertices as intersection of 2 edges
// static int Tve[4][2] = {{0,1},{0,3},{1,3},{2,4}};
static int Tve[4][2] = {{0,2},{0,1},{1,2},{3,4}};

// gives edges as intersection of 2 faces
// static int Tef[6][2] = {{0,1},{0,2},{1,2},{0,3},{1,3},{2,3}};
static int Tef[6][2] = {{0,1},{0,2},{0,3},{1,3},{1,2},{2,3}};

// gives vertices as intersection of 3 faces
// static int Tvf[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
static int Tvf[4][3] = {{0,1,3},{0,1,2},{0,2,3},{1,2,3}};

mEntity *mTet::getTemplate(int ith, int what, int with)const
{
  switch(what)
    {
    case 0:
      if(theAdjacencies[0] && theAdjacencies[0]->size()==4 )return (*theAdjacencies[0])[ith];
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
      if(theAdjacencies[2])
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
      if(with == 0)
	{
	  if (isAdjacencyCreated (0))
	    {
	      return new mFace ((mVertex*)theAdjacencies[0]->get(Tfv[ith][0]),
				(mVertex*)theAdjacencies[0]->get(Tfv[ith][1]),
				(mVertex*)theAdjacencies[0]->get(Tfv[ith][2]),
				theClassification);
	    }
	  else
	    {
	      return new mFace ((mVertex*)get(0,Tfv[ith][0]),
				(mVertex*)get(0,Tfv[ith][1]),
				(mVertex*)get(0,Tfv[ith][2]),
				theClassification);
	    }
	}

      else if (with == 1) 
	return new mFace ((mEdge*)get(1,Tfe[ith][0]),
			  (mEdge*)get(1,Tfe[ith][1]),
			  (mEdge*)get(1,Tfe[ith][2]),
			  theClassification);
    default : throw new mException (__LINE__,__FILE__,"this tet template is not done");
    }
  return nullptr;
}

mTet::~mTet()
= default;

} // end of namespace
