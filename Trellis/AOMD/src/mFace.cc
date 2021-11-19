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
#include "mFace.h"
#include "mEntityContainer.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mUnits.h"
#include "mException.h"

#include <cstdio>

using std::set;
using std::cout;

namespace AOMD {

mFace::mFace()
  : mEntity()
{
}

mFace::mFace(mVertex *v1, mVertex *v2, mVertex *v3,GEntity *classification)
  : mEntity ()
{
  typ = TRI;
  theAdjacencies[0] = new mAdjacencyContainer (3);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
  theClassification = classification;
  iD = v1->getRAND()+ v2->getRAND() + v3->getRAND();
}

mFace::mFace(mEdge *e1, mEdge *e2, mEdge *e3,GEntity *classification, int *dir)
  : mEntity ()
{
  typ = TRI;
  theAdjacencies[1] = new mAdjacencyContainer (3);
  if (!dir)
    {
      theAdjacencies[1]->add(e1);
      theAdjacencies[1]->add(e2);	
      theAdjacencies[1]->add(e3);	
    }
  else
    {
      theAdjacencies[1]->add(e1);
      mVertex *v0=(mVertex*)e1->get(0,0);
      if( dir[0]>0 ) 
	{
	  if( (mVertex*)e2->get(0,0)==v0 || (mVertex*)e2->get(0,1)==v0 )
	    {
	      theAdjacencies[1]->add(e3);
	      theAdjacencies[1]->add(e2);
	    }
	  else
	    {
	      theAdjacencies[1]->add(e2);
	      theAdjacencies[1]->add(e3);
	    }
	}
      else
	{
	  if( (mVertex*)e2->get(0,0)==v0 || (mVertex*)e2->get(0,1)==v0 )
	    {
	      theAdjacencies[1]->add(e2);
	      theAdjacencies[1]->add(e3);
	    }
	  else
	    {
	      theAdjacencies[1]->add(e3);
	      theAdjacencies[1]->add(e2);
	    }
	}
    }

  theClassification = classification;

  mVertex *v1 = (mVertex*)get(0,0);
  mVertex *v2 = (mVertex*)get(0,1);
  mVertex *v3 = (mVertex*)get(0,2);
  if(!v1 || !v2 || !v3)
    {
      char text[256];
      sprintf(text,"error in the creation of a triangle with its edges\n edges (%d %d), (%d %d) and (%d %d) cannot be nested",
	      e1->vertex(0)->getId(),e1->vertex(1)->getId(),
	      e2->vertex(0)->getId(),e2->vertex(1)->getId(),
	      e3->vertex(0)->getId(),e3->vertex(1)->getId());
      throw new mException (__LINE__,__FILE__,text);
    }
  iD = v1->getRAND() + v2->getRAND() + v3->getRAND(); 
#ifndef FLEXDB
  theAdjacencies[0] = new mAdjacencyContainer (3);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
#endif  
}

mFace::mFace(mVertex *v1, mVertex *v2, mVertex *v3, mVertex *v4, GEntity *classification)
  : mEntity ()
{
  typ = QUAD;
  theAdjacencies[0] = new mAdjacencyContainer (4);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
  theAdjacencies[0]->add(v4);	
  theClassification = classification;
  iD = v1->getRAND()+ v2->getRAND() + v3->getRAND() + v4->getRAND();
}

mFace::mFace(mEdge *e1, mEdge *e2, mEdge *e3, mEdge *e4, GEntity *classification, int *dir)
   : mEntity ()
{
  typ = QUAD; 

  theAdjacencies[1] = new mAdjacencyContainer (4);

  if (!dir || dir[0] > 0)
    {
      theAdjacencies[1]->add(e1);
      theAdjacencies[1]->add(e2);	
      theAdjacencies[1]->add(e3);	
      theAdjacencies[1]->add(e4);	
    }
  else
    {
      theAdjacencies[1]->add(e4);
      theAdjacencies[1]->add(e3);	
      theAdjacencies[1]->add(e2);	
      theAdjacencies[1]->add(e1);	
    }


  theClassification = classification;
  /*this will give the double of the result*/
  iD = e1->getId() + e2->getId() + e3->getId() + e4->getId(); 
  iD /=2;

  mVertex *v1 = (mVertex*)get(0,0);
  mVertex *v2 = (mVertex*)get(0,1);
  mVertex *v3 = (mVertex*)get(0,2);
  mVertex *v4 = (mVertex*)get(0,3);
  if(!v1 || !v2 || !v3 || !v4)
    {
      char text[256];
      sprintf(text,"\terror in the creation of a quad with its edges\n\tedges (%d %d), (%d %d), (%d %d) and (%d %d) cannot be nested",
	      e1->vertex(0)->getId(),e1->vertex(1)->getId(),
	      e2->vertex(0)->getId(),e2->vertex(1)->getId(),
	      e3->vertex(0)->getId(),e3->vertex(1)->getId(),
	      e4->vertex(0)->getId(),e4->vertex(1)->getId()
	      );
      throw new mException (__LINE__,__FILE__,text);
    }
}

int mFace:: getNbTemplates (int what) const
{
  return(typ == TRI)?3:4;
}

mEntity *mFace::getTemplate(int ith, int what, int with)const
{
  const int nbT = (typ == TRI)?3:4;

  switch(what)
    {
    case 0:
      if(theAdjacencies[0])return (*theAdjacencies[0])[ith];
      if(theAdjacencies[1])
      {
          mEdge *e1 = (mEdge*)get(1,(ith+nbT-1) % nbT);
          mEdge *e2 = (mEdge*)get(1,ith);
          return e1->commonVertex(e2);
      }
      break;
    case 1:
      if(isAdjacencyCreated(0))
      {
         return new mEdge ((mVertex*)theAdjacencies[0]->get(ith),
			    (mVertex*)theAdjacencies[0]->get((ith+1)%(nbT)),
			    theClassification);
      }
      else
      {
         return new mEdge ((mVertex*)get(0,ith),
			    (mVertex*)get(0,(ith+1)%(nbT)),
			    theClassification);
      }
      break;
    default:
      cout << "mFace::getTemplate error : level " << what << " is incorrect\n";
      throw 1;
      return nullptr;
    }
  return nullptr;
}

mVertex * mFace::commonVertex (mFace *f1, mFace *f2)
{
  mVertex * thisVertices[5];
  mVertex * f1Vertices[5];
  mVertex * f2Vertices [5];

  const int thisSize = getNbTemplates(0);
  const int f1Size   = f1->getNbTemplates(0);
  const int f2Size   = f2->getNbTemplates(0);

{  for(int i=0;i<thisSize;i++)thisVertices[i] = (mVertex*)get(0,i);}
{  for(int i=0;i<f1Size;i++)f1Vertices[i] = (mVertex*)f1->get(0,i);}
{  for(int i=0;i<f2Size;i++)f2Vertices[i] = (mVertex*)f2->get(0,i);}
    
  for( int i=0;i<thisSize;i++)
    {
      if(std::find(f1Vertices,f1Vertices+f1Size,thisVertices[i]) != f1Vertices+f1Size &&
	 std::find(f2Vertices,f2Vertices+f2Size,thisVertices[i]) != f2Vertices+f2Size)
	return thisVertices[i];
    }
  return nullptr;
}

mEdge * mFace::commonEdge (mFace *f1)
{
  if(!isAdjacencyCreated(1))return nullptr;

  for(int i=0;i<size(1);i++)
    {
      mEntity *e1 = get(1,i);
      for(int j=0;j<f1->size(1);j++)
	{
	  if(e1 == f1->get(1,j))return (mEdge*)e1;
	}      
    }
  return nullptr;
}

mFace::~mFace()
= default;

int mFace::getLevel()const
{
  return 2;
}

} // end of namespace
