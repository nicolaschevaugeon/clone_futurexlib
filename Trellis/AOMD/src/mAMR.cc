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
#include <cstdio>
#include "mMesh.h"
#include "mVertex.h"
#include "mMirrorEntity.h"
#include "mFace.h"
#include "mEdge.h"
#include "mHex.h"
#include "mTet.h"
#include "mAOMD.h"
#include "ParUtil.h"
#include "AOMD_OwnerManager.h"
#include "mExchangeData.h"
#include "mBuildAdj.h"

#include <iostream>
#include <list>
#include <set>

using std::cout;
using std::list;
using std::set;



extern "C" {
void MeshModel_middlePoint (AOMD::mMesh* m, AOMD::mEntity *e, double *p);
}

namespace AOMD {

//*****************************************
// functions moved from gEntity
  Trellis_Util::mPoint mMesh::COG (mEntity *e) 
  {
    Trellis_Util::mPoint p(0,0,0);
    /// we have a model ...
#ifndef SIM
    if (e->getLevel() == 1) 
      {
	double pd[3];
	MeshModel_middlePoint (this,e,pd);
	return Trellis_Util::mPoint (pd[0],pd[1],pd[2]);
      }
    /// centriod ...
    else
#endif
      {
	for(int i=0;i<e->getNbTemplates(0);i++)
	  {
	    p += (((mVertex*)e->get(0,i))->point());
	  }
	return (p*(1./(double)e->getNbTemplates(0)));	
      }
  }


static void  createDownward(mMesh *m, mEntity *e , int i, int with)
{
  createDownwardFunctor f(i,with,m);
  f(e);
}

template <class T>
void myswap (T &t1, T &t2)
{
  T temp;
  temp = t1;
  t1 = t2;
  t2 = temp;
}

int mMesh::getRefinementLevel(mEntity *e)
{
  if(!e->parent())return 0;
  else return 1 + getRefinementLevel(e->parent());
}


int mMesh::getRefinementDepth(mEntity *e)
{
  int dim = e->getLevel();
  if(e->isAdjacencyCreated(dim))
    {
      int theMax=0;
      for(int i=0;i<e->size(dim);i++)
	{
	  mEntity *sub = e->get(dim,i);
	  int x  =  getRefinementDepth(sub) + 1;
	  theMax = (theMax>x)?theMax:x;
	}
      
       return theMax;
    }
  else return 0;
}

static mVertex * splitEdge(mMesh *m, mEdge *e, mEntity *src)
{

  if(e->isAdjacencyCreated(1))return ((mEdge*)e->get(1,0))->commonVertex((mEdge*)e->get(1,1));

  mVertex *newv;
  GEntity *g = e->getClassification();
  Trellis_Util::mPoint pp (m->COG(e));

  if(e->getData(AOMD_Util::Instance()->getAtt1()))
    {  
      //      printf("proc %d edge splitted with vertex %d ",
      //	      ParUtil::Instance()->rank(), e->getAttachedInt(AOMD_Util::Instance()->getAtt1()));e->print();
      newv = m->createVertex (e->getAttachedInt(AOMD_Util::Instance()->getAtt1()),pp(0),pp(1),pp(2),g);
    }
  else
    {
      e->print();
      ParUtil::Instance()->Msg(ParUtil::ERROR,"Unflagged edge to split ");
      newv = m->createVertex (pp(0),pp(1),pp(2),g);
    }

  if(e->getData(AOMD_Util::Instance()->getAtt2()))
    {
      //      printf("creating new mirror vertex %d \n",e->getAttachedInt(AOMD_Util::Instance()->getAtt2()));
      mMirrorVertex *mv =  (mMirrorVertex*)m->getVertex(e->getAttachedInt(AOMD_Util::Instance()->getAtt2()));
      //jf a faire if(!mv)mv = m->createMirrorVertex(e->getAttachedInt(AOMD_Util::Instance()->getAtt2()));
      //jf a faire if(mv->getType() != mEntity::MIRROR)
      //jf a faire ParUtil::Instance()->Msg(ParUtil::ERROR,"Not a mirror vertex\n");
      mv->addCopy(newv);
      //      mv->print();
    }

  mEdge *e1 = m->createEdge(e->vertex(0),newv,e->getClassification());
  mEdge *e2 = m->createEdge(newv,e->vertex(1),e->getClassification());
  
  e1->setParent(e);
  e2->setParent(e);

  e->add(e1);e->add(e2);  

  e1->setCommonBdry(e->getCommonBdry());
  e2->setCommonBdry (e->getCommonBdry());

  return newv;
}


mVertex *mMesh::splitQuad(mFace *e, mEntity *src)
{
  mMesh *m = this;

  if(e->isAdjacencyCreated(2))return ((mFace*)e->get(2,0))->commonVertex
				((mFace*)e->get(2,1),(mFace*)e->get(2,2));
  
  GEntity *g = e->getClassification();
  
  mEdge *e1 = (mEdge*)e->get(1,0);
  mEdge *e2 = (mEdge*)e->get(1,1);
  mEdge *e3 = (mEdge*)e->get(1,2);
  mEdge *e4 = (mEdge*)e->get(1,3);
  splitEdge(m,e1,src);
  splitEdge(m,e2,src);
  splitEdge(m,e3,src);
  splitEdge(m,e4,src);  
  mVertex *v2 = e1->commonVertex(e2);
  mVertex *v3 = e2->commonVertex(e3);
  mVertex *v4 = e3->commonVertex(e4);
  mVertex *v1 = e4->commonVertex(e1);
   
  mEdge *e11 = (mEdge*)e1->get(1,0);
  mEdge *e12 = (mEdge*)e1->get(1,1);

  if(!e11->commonVertex(e4))
    {
      myswap<mEdge*> (e11,e12);
    }

  mEdge *e21 = (mEdge*)e2->get(1,0);
  mEdge *e22 = (mEdge*)e2->get(1,1);

  if(!e21->commonVertex(e1))
    {
      myswap<mEdge*> (e21,e22);
    }

  mEdge *e31 = (mEdge*)e3->get(1,0);
  mEdge *e32 = (mEdge*)e3->get(1,1);

  if(!e31->commonVertex(e2))
    {
      myswap<mEdge*> (e31,e32);
    }

  mEdge *e41 = (mEdge*)e4->get(1,0);
  mEdge *e42 = (mEdge*)e4->get(1,1);

  if(!e41->commonVertex(e3))
    {
      myswap<mEdge*> (e41,e42);
    }

  mVertex *v12 = e11->commonVertex(e12);
  mVertex *v23 = e21->commonVertex(e22);
  mVertex *v34 = e31->commonVertex(e32);
  mVertex *v41 = e41->commonVertex(e42);

  Trellis_Util::mPoint pp (m->COG(e));

  mVertex *vnew;
  
  if(e->getData(AOMD_Util::Instance()->getAtt1()))
    {  
      /*      printf("proc %d quad splitted with vertex %d ",
	      ParUtil::Instance()->rank(), e->getAttachedInt(AOMD_Util::Instance()->getAtt1()));e->print();*/
      vnew = createVertex (e->getAttachedInt(AOMD_Util::Instance()->getAtt1()),pp(0),pp(1),pp(2),g);
    }
  else
    {
      ParUtil::Instance()->Msg(ParUtil::WARNING,"Unflagged quad to split\n");
      vnew = createVertex (pp(0),pp(1),pp(2),e->getClassification());  
    }
  
  if(e->getData(AOMD_Util::Instance()->getAtt2()))
    {
      //      printf("creating new mirror vertex %d \n",e->getAttachedInt(AOMD_Util::Instance()->getAtt2()));
      mMirrorVertex *mv =  (mMirrorVertex*)getVertex(e->getAttachedInt(AOMD_Util::Instance()->getAtt2()));
      //jf a faire if(!mv)mv = createMirrorVertex(e->getAttachedInt(AOMD_Util::Instance()->getAtt2()));
      mv->addCopy(vnew);
    }


  mFace *f1 = createFaceWithVertices(v1,v12,vnew,v41,e->getClassification());
  mFace *f2 = createFaceWithVertices(v2,v23,vnew,v12,e->getClassification());
  mFace *f3 = createFaceWithVertices(v3,v34,vnew,v23,e->getClassification());
  mFace *f4 = createFaceWithVertices(v4,v41,vnew,v34,e->getClassification());

  f1->setParent(e);
  f2->setParent(e);
  f3->setParent(e);
  f4->setParent(e);

  createDownward(this,f1,1,0);
  createDownward(this,f2,1,0);
  createDownward(this,f3,1,0);
  createDownward(this,f4,1,0);
  
  e->add(f1);
  e->add(f2);
  e->add(f3);
  e->add(f4);
  
  if(!find(v12))add(v12);
  if(!find(v23))add(v23);
  if(!find(v34))add(v34);
  if(!find(v41))add(v41);

  f1->setCommonBdry (e->getCommonBdry());
  f2->setCommonBdry (e->getCommonBdry());
  f3->setCommonBdry (e->getCommonBdry());
  f4->setCommonBdry (e->getCommonBdry());


  return vnew;
}

void mMesh::splitTriangle(mFace *e, mEntity *src)
{
  mMesh *m = this;
  if(e->isAdjacencyCreated(2))return;

  mEdge *e1 = (mEdge*)e->get(1,0);
  mEdge *e2 = (mEdge*)e->get(1,1);
  mEdge *e3 = (mEdge*)e->get(1,2);
  splitEdge(m,e1,src);
  splitEdge(m,e2,src);
  splitEdge(m,e3,src);
  mVertex *v2 = e1->commonVertex(e2);
  mVertex *v3 = e2->commonVertex(e3);
  mVertex *v1 = e3->commonVertex(e1);
   
  mEdge *e11 = (mEdge*)e1->get(1,0);
  mEdge *e12 = (mEdge*)e1->get(1,1);

  if(!e11->commonVertex(e3))
    {
      myswap<mEdge*> (e11,e12);
    }

  mEdge *e21 = (mEdge*)e2->get(1,0);
  mEdge *e22 = (mEdge*)e2->get(1,1);

  if(!e21->commonVertex(e1))
    {
      myswap<mEdge*> (e21,e22);
    }

  mEdge *e31 = (mEdge*)e3->get(1,0);
  mEdge *e32 = (mEdge*)e3->get(1,1);

  if(!e31->commonVertex(e2))
    {
      myswap<mEdge*> (e31,e32);
    }

  mVertex *v12 = e11->commonVertex(e12);
  mVertex *v23 = e21->commonVertex(e22);
  mVertex *v31 = e31->commonVertex(e32);
  
  mEdge *enew1 = createEdge (v12,v23,e->getClassification());
  mEdge *enew2 = createEdge (v23,v31,e->getClassification());
  mEdge *enew3 = createEdge (v31,v12,e->getClassification());
  
  mFace *f1 = createFaceWithVertices (v1,v12,v31,e->getClassification());
  mFace *f2 = createFaceWithVertices (v12,v2,v23,e->getClassification());
  mFace *f3 = createFaceWithVertices (v31,v23,v3,e->getClassification());
  mFace *f4 = createFaceWithVertices (v12,v23,v31,e->getClassification());

  f1->setParent(e);
  f2->setParent(e);
  f3->setParent(e);
  f4->setParent(e);
  //  int r =  getRefinementLevel(e) + 1;
  //  attachRefinementLevel( f1, r);
  //  attachRefinementLevel( f2, r);
  //  attachRefinementLevel( f3, r);
  //  attachRefinementLevel( f4, r);

  f1->add(e11)  ;  f1->add(enew3); f1->add(e32);
  f2->add(e12)  ;  f2->add(e21)  ; f2->add(enew1);
  f3->add(enew2);  f3->add(e22)  ; f3->add(e31);
  f4->add(enew1);  f4->add(enew2); f4->add(enew3);
    
  e->add(f1);
  e->add(f2);
  e->add(f3);
  e->add(f4);


  f1->setCommonBdry(e->getCommonBdry());
  f2->setCommonBdry(e->getCommonBdry());
  f3->setCommonBdry(e->getCommonBdry());
  f4->setCommonBdry(e->getCommonBdry());

}

void mMesh::unsplitTriangle(mFace *e)
{
  mFace *f1 = (mFace*)e->get(2,0);
  mFace *f2 = (mFace*)e->get(2,1);
  mFace *f3 = (mFace*)e->get(2,2);
  mFace *f4 = (mFace*)e->get(2,3);

  f1->deleteParent();
  f2->deleteParent();
  f3->deleteParent();
  f4->deleteParent();

  mEdge *e1 = nullptr,*e2 = nullptr,*e3 = nullptr;

  if(!f1->commonEdge(f3) && !f1->commonEdge(f2))
    {
      e1 = f4->commonEdge(f1);
      e2 = f4->commonEdge(f2);
      e3 = f4->commonEdge(f3);
    }
  else if(!f1->commonEdge(f2) && !f1->commonEdge(f4))
    {
      e1 = f3->commonEdge(f1);
      e2 = f3->commonEdge(f2);
      e3 = f3->commonEdge(f4);
    }
  else if(!f1->commonEdge(f4) && !f1->commonEdge(f3))
    {
      e1 = f2->commonEdge(f1);
      e2 = f2->commonEdge(f3);
      e3 = f2->commonEdge(f4);
    }
  else if(!f2->commonEdge(f3) && !f2->commonEdge(f4))
    {
      e1 = f1->commonEdge(f2);
      e2 = f1->commonEdge(f3);
      e3 = f1->commonEdge(f4);
    }

  if(!e1 || !e2 || !e3)
    {
      printf("error in unsplitTre, %p %p %p\n",e1,e2,e3);  
      printf("f1 %p %p %p\n",f1->get(1,0),f1->get(1,1),f1->get(1,2));
      printf("f2 %p %p %p\n",f2->get(1,0),f2->get(1,1),f2->get(1,2));
      printf("f3 %p %p %p\n",f3->get(1,0),f3->get(1,1),f3->get(1,2));
      printf("f4 %p %p %p\n",f4->get(1,0),f4->get(1,1),f4->get(1,2));
      throw 1;
    }

  DEL (e1);
  DEL (e2);
  DEL (e3);

  DEL(f1);
  DEL(f2);
  DEL(f3);
  DEL(f4);

  e->deleteAdjacencies(2);
}

void mMesh::unsplitQuad(mFace *e)
{
  mFace *f1 = (mFace*)e->get(2,0);
  mFace *f2 = (mFace*)e->get(2,1);
  mFace *f3 = (mFace*)e->get(2,2);
  mFace *f4 = (mFace*)e->get(2,3);
  f1->deleteParent();
  f2->deleteParent();
  f3->deleteParent();
  f4->deleteParent();

  mEdge *e1 = nullptr,*e2 = nullptr,*e3 = nullptr,*e4 = nullptr;

  if(!f1->commonEdge(f3))
    {
      e1 = f1->commonEdge(f2);
      e2 = f1->commonEdge(f4);
      e3 = f3->commonEdge(f2);
      e4 = f3->commonEdge(f4);
    }
  else if (!f1->commonEdge(f4))
    {
      e1 = f1->commonEdge(f2);
      e2 = f1->commonEdge(f3);
      e3 = f4->commonEdge(f2);
      e4 = f4->commonEdge(f3);
    }
  else if (!f1->commonEdge(f2))
    {
      e1 = f1->commonEdge(f3);
      e2 = f1->commonEdge(f4);
      e3 = f2->commonEdge(f3);
      e4 = f2->commonEdge(f4);
    }
  if(!e1 || !e2 || !e3 || !e4)
    {
      printf("error in unsplitQuad, %p %p %p %p\n",e1,e2,e3,e4);  
      printf("f1 %p %p %p %p\n",f1->get(1,0),f1->get(1,1),f1->get(1,2),f1->get(1,3));
      printf("f2 %p %p %p %p\n",f2->get(1,0),f2->get(1,1),f2->get(1,2),f2->get(1,3));
      printf("f3 %p %p %p %p\n",f3->get(1,0),f3->get(1,1),f3->get(1,2),f3->get(1,3));
      printf("f4 %p %p %p %p\n",f4->get(1,0),f4->get(1,1),f4->get(1,2),f4->get(1,3));
    }

  DEL (e1);
  DEL (e2);
  DEL (e3);
  DEL (e4);
  mVertex *v = f4->commonVertex(f1,f2);

  DEL(f1);
  DEL(f2);
  DEL(f3);
  DEL(f4);
  DEL(v);

  e->deleteAdjacencies(2);
}

void mMesh::unsplitHex(mHex *h)
{
  mHex *H[9];
  mFace *F[13]; // because we evaluate F[12] !! bad mistake that I've done, forgive me...
  mEdge *E[7];
  mVertex *V;
  int i;

  //  printf("unsplitting HEX ");h->print();
  
  mAdjacencyContainer c;

  if(!h->isAdjacencyCreated(3))return;
  
  
  for(i=0;i<8;i++)
    {
      H[i] = (mHex*)h->get(3,i);      
      //      printf("sub Hex(%d) = ",i);H[i]->print();
    }

  int NF = 0,NE=0;
  for(i=0;i<8;i++)
    {
      for(int j=0;j<i;j++)
	{
	  if((F[NF] = H[i]->commonFace(H[j])))NF++;

	}
    }

  for(i=0;i<NF;i++)
    {
      for(int j=0;j<i;j++)
	{
	  //	  F[i]->print();F[j]->print();
	  if((E[0]=F[i]->commonEdge(F[j])))
	    {
	      //	      E[0]->print();
	      c.appendUnique(E[0]);
	    }
	}
    }
  


  //  printf("%d faces to delete\n",NF);

  for(i=0;i<8;i++)DEL(H[i]);
  for(i=0;i<NF;i++)
    {
      //      printf("Internal Face(%d) = ",i);F[i]->print();
      DEL(F[i]);
    }
  NE = c.size();
  //  printf("%d edges to delete\n",NE);
{  for(int i=0;i<NE;i++)
    {
      E[i] = (mEdge*)c[i];
    }
}
  V = E[0]->commonVertex(E[1]);

  for(i=0;i<NE;i++)
    {
      //      printf("Internal Edge(%d) = ",i);E[i]->print();
      DEL(E[i]);
    }

  //  printf("internal Vertex = ");V->print();
  DEL(V);
  h->deleteAdjacencies(3);
}

void mMesh::unsplitTet(mTet *t)
{
  mTet *T[8];
  mTet *TC[4];
  mFace *F[23]; 
  mEdge *theE;
  int i,k=0;
  
  if(!t->isAdjacencyCreated(3))return;
  
  for(i=0;i<8;i++)
    {
      T[i] = (mTet*)t->get(3,i);      
      if(!T[i]->find(t->get(0,0)) && 
	 !T[i]->find(t->get(0,1)) && 
	 !T[i]->find(t->get(0,2)) && 
	 !T[i]->find(t->get(0,3)))TC[k++] = T[i];
    }
  
  if(k!=4)printf("another major error connor\n");

  int NF = 0;
  for(i=0;i<8;i++)
    {
      for(int j=0;j<i;j++)
	{
	  if((F[NF] = T[i]->commonFace(T[j])))NF++;
	}
    }  

  theE = TC[0]->commonEdge (TC[1],TC[2]);

  if(!theE)
    {
      printf("Error major connor\n");
    }
  for(i=0;i<8;i++)DEL(T[i]);
  for(i=0;i<NF;i++)
    {
      DEL(F[i]);
    }
  DEL(theE);
  t->deleteAdjacencies(3);
}



mVertex * mMesh::splitHex(mHex *h)
{
  mMesh *m = this;
  if(h->isAdjacencyCreated(3))
    {
      printf("there is something wrong in split hex !!!!\n");
      return nullptr; // wrong !!
    }
  int i;
  GEntity *g = h->getClassification();
  mFace *f[6];
  mEdge *e[12];
  mVertex *v[8];
  mVertex *vf[6];
  mVertex *ve[12];
  mHex *H[8];

  //  printf("splitting an hex ");h->print();

  mVertex *vnew;
  Trellis_Util::mPoint pp(m->COG(h));
  //  printf("cog : %f %f %f\n",pp(0),pp(1),pp(2));
  if(h->getData(AOMD_Util::Instance()->getAtt1()))
    {  
      
      /*      printf("proc %d hex splitted with vertex %d ",
	      ParUtil::Instance()->rank(), h->getAttachedInt(AOMD_Util::Instance()->getAtt1()));h->print();*/
      vnew = createVertex (h->getAttachedInt(AOMD_Util::Instance()->getAtt1()),pp(0),pp(1),pp(2),g);
    }
  else
    vnew = createVertex (pp(0),pp(1),pp(2),g);
  
  for(i=0;i<8;i++)
    {
      v[i] = (mVertex*)h->get(0,i);
    }
  for(i=0;i<12;i++)
    {
      e[i] = (mEdge*)h->get(1,i);
      //      printf("splitting edge ");e[i]->print();
      ve[i] = splitEdge(m,e[i],h);
    }
  for(i=0;i<6;i++)
    {
      f[i] = (mFace*)h->get(2,i);
      //      printf("splitting face ");f[i]->print();
      vf[i] = splitQuad(f[i],h);
    }
  
  H[0] = createHexWithVertices ( v[0],ve[0],vf[0],ve[3],
				 ve[4],vf[1],vnew,vf[4],g);
  H[1] = createHexWithVertices ( ve[0],v[1],ve[1],vf[0],
				 vf[1],ve[5],vf[2],vnew,g);
  H[2] = createHexWithVertices ( ve[4],vf[1],vnew,vf[4],
				 v[4],ve[8],vf[5],ve[11],g);
  H[3] = createHexWithVertices ( vf[1],ve[5],vf[2],vnew,
				     ve[8],v[5],ve[9],vf[5],g);
				     
  H[4] = createHexWithVertices ( ve[3],vf[0],ve[2],v[3],
				 vf[4],vnew,vf[3],ve[7],g);
  H[5] = createHexWithVertices ( vf[0],ve[1],v[2],ve[2],
				 vnew,vf[2],ve[6],vf[3],g);
  H[6] = createHexWithVertices ( vf[4],vnew,vf[3],ve[7],
				 ve[11],vf[5],ve[10],v[7],g);
  H[7] = createHexWithVertices ( vnew,vf[2],ve[6],vf[3],
				 vf[5],ve[9],v[6],ve[10],g);
				     
  //  printf("hex splitted !\n");
  for(i=0;i<8;i++)
    {
      // if we do not want r->e adjacencices stored
      if(h->isAdjacencyCreated(1))createDownward(this,H[i],1,0);
      // create faces with vertices
      createDownward(this,H[i],2,0);
      // created edges of faces
      for(int j=0;j<6;j++)
      {
	createDownward(this,H[i]->get(2,j),1,0);
      }
      h->add(H[i]);
      H[i]->setParent(h);
      if(!h->isAdjacencyCreated(0))H[i]->deleteAdjacencies(0);
    }

  return vnew;
}



void mMesh::splitTet(mTet *t)
{
  //  printf("proc %d split tet ",ParUtil::Instance()->rank()); t->print();
  mMesh *m = this;
  if(t->isAdjacencyCreated(3))
    {
      printf("there is something wrong in split tet !!!!\n");
      return;
    }
  int i;
  GEntity *g = t->getClassification();
  mFace *f[4];
  mEdge *e[6];
  mVertex *v[4];
  mVertex *ve[6];
  mTet *T[8];
  int edge_permutation[6] = {0,2,3,1,4,5};

  for(i=0;i<4;i++)
    {
      v[i] = (mVertex*)t->get(0,i);
    }
  for(i=0;i<6;i++)
    {
      e[i] = (mEdge*)t->get(1,edge_permutation[i]);
      ve[i] = splitEdge(m,e[i],t);
    }
  for(i=0;i<4;i++)
    {
      f[i] = (mFace*)t->get(2,i);
      splitTriangle(f[i],t);
    }

  // four corners  

  T[0] = createTetWithVertices ( v[0],ve[0],ve[2],ve[1],g);
  T[1] = createTetWithVertices ( v[1],ve[0],ve[3],ve[4],g);
  T[2] = createTetWithVertices ( v[2],ve[3],ve[1],ve[5],g);
  T[3] = createTetWithVertices ( v[3],ve[2],ve[4],ve[5],g);
				     
  // four inside (internal edge is chosen ve[3] -> ve[2])
  T[4] = createTetWithVertices ( ve[3],ve[5],ve[2],ve[4],g);
  T[5] = createTetWithVertices ( ve[3],ve[2],ve[0],ve[4],g);
  T[6] = createTetWithVertices ( ve[2],ve[5],ve[3],ve[1],g);
  T[7] = createTetWithVertices ( ve[0],ve[2],ve[3],ve[1],g);

  for(i=0;i<8;i++)
    {
      // if we do not want r->e adjacencices stored
      if(t->isAdjacencyCreated(1))createDownward(this,T[i],1,0);
      // create faces with vertices
      createDownward(this,T[i],2,0);
      // created edges of faces
      for(int j=0;j<4;j++)
      {
	createDownward(this,T[i]->get(2,j),1,0);
      }
      t->add(T[i]);
      T[i]->setParent(t);
      if(!t->isAdjacencyCreated(0))T[i]->deleteAdjacencies(0);
    }
}


static void unsplitEdge(mMesh *m, mEdge *e)
{
  mEdge *e1 = (mEdge*)e->get(1,0);
  mEdge *e2 = (mEdge*)e->get(1,1);
  e1->deleteParent();
  e2->deleteParent();
  mVertex *v = e1->commonVertex(e2);
  m->DEL(e1);
  m->DEL(e2);
  m->DEL(v);
  e->deleteAdjacencies(1);
}

void mMesh::unsplit(list<mEntity*> &toUnsplit, mSplitCallbacks &f)
{
   
  //  printf("phase 1...\n");

  f.unsplitCallback(toUnsplit);

  for(std::list<mEntity*>::const_iterator it = toUnsplit.begin();
      it != toUnsplit.end();
      ++it)
    {
      mEntity *e = *it;
      switch(e->getType())
	{
	case mEntity::VERTEX:
	  break;
	case mEntity::EDGE:
	  unsplitEdge(this,(mEdge*)e);
	  break;
	case mEntity::TRI:
	  unsplitTriangle((mFace*)e);
	  break;      
	case mEntity::QUAD:
	  unsplitQuad((mFace*)e);
	  break;      
	case mEntity::HEX:
	  unsplitHex((mHex*)e);
	  break;      
	case mEntity::TET:
	  unsplitTet((mTet*)e);
	  break;      
	default:
	  std::cout<< "Error in " << __FILE__<< ":" << __LINE__<<  " case " << e->getType() << " not treated in switch" << std::endl;
	  throw;
	}
    }
}

inline void markEntity (mEntity *e , bool mark)
{
  if(mark)
    {
      if(!e->getData(AOMD_Util::Instance()->getAtt3()))
	e->attachInt(AOMD_Util::Instance()->getAtt3(),1);
    }
  else e->deleteData(AOMD_Util::Instance()->getAtt3());
}

/*
  If an entity is used on one processor, we mark it
  used for all remote copies.
*/

class Used_DataExchanger : public AOMD_DataExchanger
{
  int dim;
public :
  Used_DataExchanger(int d):dim(d){}
  int tag() const override {return 184;}
  void * AP_alloc_and_fill_buffer (mEntity *e, AOMD_SharedInfo &si, int) override;
  void receiveData (int pid, void *buf) override;   
};

void * Used_DataExchanger :: AP_alloc_and_fill_buffer (mEntity *e, 
						       AOMD_SharedInfo &si, int d_tag)
{
  if(e->getLevel() != dim)return nullptr;
  if(e->getData(AOMD_Util::Instance()->getAtt3()))
    {

      void *buf=malloc(sizeof(rp_int));

      rp_int *rp = (rp_int*)buf;
      rp->entity =  si.getRemotePointer();
      rp->i = 1; 
      return buf;
    }
  else return nullptr;
}


void Used_DataExchanger :: receiveData (int from, void *buf)
{
  rp_int *castbuf= (rp_int *)buf;  
  markEntity(castbuf->entity,true);
}

/*-----------------------------------------------------------------*/

void markUsedEntities (mMesh *theMesh, int dim, bool mark, bool all = false)
{

  //  printf("size(%d) = %d\n",dim,size(dim));

  for(mMesh::iterall it = theMesh->beginall(dim);
      it != theMesh->endall(dim);
      ++it)
    {
      mEntity *ent = *it;
      
      if(all)
	{
	  set<mEntity*> theWholeFamily;
	  ent->getAllEntitiesOfAllDimensions(theWholeFamily);
	  for(std::set<mEntity*>::iterator it2 = 
		theWholeFamily.begin(); it2 != theWholeFamily.end();++it2)
	    {
	      markEntity(*it2,mark);
	    }	  
	}
      else
	{      
	  for(int j=0;j<dim;j++)
	    {
	      if(ent->isAdjacencyCreated(j))
		{
		  for(int i=0;i<ent->size(j);i++)
		    {
		      markEntity(ent->get(j,i),mark);
		    }
		}
	    }
	}
    }
}

void mMesh::cleanup (int dim,  mSplitCallbacks &f)
{ 

  markUsedEntities(this,dim,true);

  //  if(ParUtil::Instance()->master())printf("--- cleanup %d ---\n",dim);

  while(1)
    {      
      //#ifdef PARALLEL
      Used_DataExchanger uda (dim-1);
      exchangeDataOnPartBdrys (uda);
      //#endif
      list<mEntity*> toUnsplit;
      for(mMesh::iterall it = beginall(dim-1);
	  it != endall(dim-1);
	  ++it)
	{
	  
	  mEntity *ent = *it;
	  if(getRefinementDepth(ent)==1)
	    {
	      bool haveToUnsplit = true;

	      // TEST

	      // if (ent->setCommonBdry()) haveToUnsplit = false;

	      for(int i=0;i<ent->size(dim-1);i++)
		{		  
		  if(ent->get(dim-1,i)->getData(AOMD_Util::Instance()->getAtt3()))
		    {	
		      haveToUnsplit = false;
		    }
		}	      
	      if(haveToUnsplit)toUnsplit.push_back(ent);
	    }
	}  
      int N = toUnsplit.size();

      if(toUnsplit.size())
	unsplit(toUnsplit, f);

      if(N) bdryLinkSetup();

      //      ParUtil::Instance()->Msg(ParUtil::WARNING,"%d %d enditities unsplitted \n",N,toUnsplit.size());

      if(!N)
	break;
    }
  markUsedEntities(this,dim,false,true);
  //  if(ParUtil::Instance()->master())printf("--- cleanup %d done ---\n",dim);
}

void mMesh :: cleanup_updateAdj(int dim)
{  
  markUsedEntities(this,dim,true,true);
  list<mEntity*> toDelete;
  for(int i=0;i<dim;i++)
    {
      for(mMesh::iterall it = beginall(i);
	  it != endall(i);
	  ++it)
	{
	  mEntity *ent = *it;
	  if(!ent->getData(AOMD_Util::Instance()->getAtt3()))toDelete.push_front(ent);
	}
    }
  markUsedEntities(this,dim,false,true);

  for(std::list<mEntity*>::const_iterator it = toDelete.begin();
      it != toDelete.end();
      ++it)
    {
      //mEntity *e = *it;
      DEL_updateAdj(*it);
    }
}


void mMesh :: cleanup (int dim)
{  
  markUsedEntities(this,dim,true,true);
  list<mEntity*> toDelete;
  for(int i=0;i<dim;i++)
    {
      for(mMesh::iterall it = beginall(i);
	  it != endall(i);
	  ++it)
	{
	  mEntity *ent = *it;
	  if(!ent->getData(AOMD_Util::Instance()->getAtt3()))toDelete.push_front(ent);
	}
    }
  markUsedEntities(this,dim,false,true);

  for(std::list<mEntity*>::const_iterator it = toDelete.begin();
      it != toDelete.end();
      ++it)
    {
      //mEntity *e = *it;
      DEL(*it);
    }
}

void markOnProcessorEntities (mMesh *theMesh, int dim, bool mark) 
{
 
  for(mMesh::iterall it = theMesh->beginall(dim);
      it != theMesh->endall(dim);
      ++it)
  { 
    mEntity *ent = *it; 
    set<mEntity*> family;
    ent->getAllEntitiesOfAllDimensions_oneLevel(family);
    for(std::set<mEntity*>::iterator it2 = 
        family.begin(); it2!=family.end();++it2)
    { 
      markEntity(*it2,mark); 
    } 
  } 
}


void mMesh ::getEntitiesToRemove(int dim,
                                 list<mEntity*>& vtToRemove, 
                                 list<mEntity*>& egToRemove,
                                 list<mEntity*>& fcToRemove)
{
  markOnProcessorEntities(this,dim,true);
//  list<mEntity*> toDelete;
  for(int i=0;i<dim;i++)
  {
    for(mMesh::iterall it = beginall(i); it != endall(i);++it)
    {
      mEntity *ent = *it;
      if(!ent->getData(AOMD_Util::Instance()->getAtt3()))
      {
        switch(ent->getLevel())
        {
          case 0: vtToRemove.push_front(ent); break;
          case 1: egToRemove.push_front(ent); break;
          case 2: fcToRemove.push_front(ent); break;
          default: break;
        }
      }
    }
  }

  markOnProcessorEntities(this,dim,false);
 
  vtToRemove.unique();
  egToRemove.unique();
  fcToRemove.unique();
}


/*
ref - unref procedure based on exterior criteria 
The actual version breaks existing upwards 
adjacencies
*/


void mSplitCallbacks :: splitCallback(list<mEntity*> &l)
{
  for(std::list<mEntity*>::const_iterator it = l.begin();it != l.end();++it)
    {
      splitCallback(*it);
    }
}

void mSplitCallbacks :: unsplitCallback(list<mEntity*> &l)
{
  for(std::list<mEntity*>::const_iterator it = l.begin();it != l.end();++it)
    {
      unsplitCallback(*it);
    }
}


void mMesh::refunref ( mSplitCallbacks &f )
{
  int dim = getDim();
  while(1)
    {
      list<mEntity*> splitList;
      for(iter it = begin(dim);
	  it != end(dim);
	  ++it)
	{
	  if(f(*it) == 1)splitList.push_back(*it);
	}
      ParUtil::Instance()->Msg(ParUtil::DEBUG2,"size(%d) = %d nb splits = %d\n",dim,size(dim),splitList.size());
      if(!split(splitList, f))break;
    }
  // unrefine;
  while(1)
    {
      list<mEntity*> unsplitList;
      for(iterall it = beginall(dim);it != endall(dim);++it)
	{
	  if(getRefinementDepth(*it) == 1) 
	    {
	      if(f(*it) == -1)unsplitList.push_back(*it);
	    }
	}
      int N = unsplitList.size();

      if(unsplitList.size())
      	unsplit(unsplitList,f);
      else if(!N)break;
      cleanup(dim,f);
      if(dim == 3)
	{
	  cleanup(2,f);
	}
    }
  //  printf("\n");
}

void mMesh::splitgen(list<mEntity*> &toSplit, mSplitCallbacks &f)
{
  int nbHex=0,nbQuad=0,nbEdg=0,nbTet=0,nbTri=0;
  for(std::list<mEntity*>::const_iterator it = toSplit.begin();
      it != toSplit.end();
      ++it)
    {
      mEntity *e = *it;
      switch(e->getType())
	{
	case mEntity::VERTEX:
	  break;
	case mEntity::EDGE:
	  splitEdge(this,(mEdge*)e,e);nbEdg++;
	  break;
	case mEntity::TRI:
	  splitTriangle((mFace*)e,e);nbTri++;
	  break;      
	case mEntity::QUAD:
	  splitQuad((mFace*)e,e);nbQuad++;
	  break;      
	case mEntity::HEX:
	  splitHex((mHex*)e);nbHex++;
	  break;      
	case mEntity::TET:
	  splitTet((mTet*)e);nbTet++;
	  break;  
	default:
	  std::cout<< "Error in " << __FILE__<< ":" << __LINE__<<  " case " << e->getType() << " not treated in switch" << std::endl;
	  throw;    
	}
    }
  f.splitCallback(toSplit);
  //  printf("proc %d (%d %d %d %d %d)\n",ParUtil::Instance()->rank(),nbEdg,nbQuad,nbHex,nbTet,nbTri);
}

static void setFlags (mEntity *e, int &n, bool b, AOMD_OwnerManager *o, bool per)
{  
  if(b)e->attachInt(AOMD_Util::Instance()->getAtt1(),++n);
  else e->deleteData(AOMD_Util::Instance()->getAtt1());
  if(b && per)e->attachInt(AOMD_Util::Instance()->getAtt2(),++n);
  else if(per)e->deleteData(AOMD_Util::Instance()->getAtt2());
}

static void flagEntity(mEntity *e, int &n, bool b, AOMD_OwnerManager *o)
{
  int i;

  if(b && e->isAdjacencyCreated(e->getLevel()) )return;
  bool per = o->isPeriodic(e);

  switch(e->getType())
    {
    case mEntity::EDGE:
      setFlags(e,n,b,o,per);
      //      printf("tagging an edge to %d %d %d %d %d ",e->getAttachedInt(AOMD_Util::Instance()->getAtt1()),per,b,n);e->print();
      break;
    case mEntity::QUAD:
    case mEntity::TRI:
      setFlags(e,n,b,o,per);
      for(i=0;i<e->size(1);i++)flagEntity(e->get(1,i),n,b,o);
      break;
    case mEntity::HEX:
      setFlags(e,n,b,o,per);
    case mEntity::TET:
      for(i=0;i<e->size(2);i++)flagEntity(e->get(2,i),n,b,o);
      break;
    default:
      std::cout<< "Error in " << __FILE__<< ":" << __LINE__<<  " case " << e->getType() << " not treated in switch" << std::endl;
      throw;
    }
}

class Refinement_DataExchanger : public AOMD_DataExchanger
{
  list<mEntity *> *toSplit;
  bool periodic;
  unsigned int ctag;

public :
  Refinement_DataExchanger(list<mEntity *> *l, bool per):toSplit(l),periodic(per)
  {  
    if(periodic)ctag = AOMD_Util::Instance()->getAtt2();
    else ctag = AOMD_Util::Instance()->getAtt1();
  }
  int tag() const override;
  void * AP_alloc_and_fill_buffer (mEntity *e, AOMD_SharedInfo &si, int) override;
  void receiveData (int pid, void *buf) override;   
};

int Refinement_DataExchanger :: tag () const
{
  return 112;
}

void * Refinement_DataExchanger :: AP_alloc_and_fill_buffer (mEntity *e, 
							     AOMD_SharedInfo &si, int d_tag)
{
  if(e->getData(ctag) && (si.isPeriodic() == periodic))
    {
      //      printf("sending "); e->print();si.getRemotePointer()->print();

      void *buf=malloc(sizeof(rp_int));

      rp_int *rp = (rp_int*)buf;
      rp->entity =  si.getRemotePointer();
      rp->i = e->getAttachedInt(ctag); 

      return buf;
    }
  else return nullptr;
}


void Refinement_DataExchanger :: receiveData (int from, void *buf)
{
  rp_int *castbuf= (rp_int *)buf;

  if(!castbuf->entity->getData(ctag))
    {
      toSplit->push_back(castbuf->entity);
      castbuf->entity->attachInt(ctag,castbuf->i);
    }
  else
    {

      int actid = castbuf->entity->getAttachedInt(ctag);

      if(castbuf->i < actid)
	{
	  castbuf->entity->attachInt(ctag,castbuf->i);
	}
    }
}

/*
  flag periodic counterparts
  if entity is to be split and is periodic, 
  we must flag all remote entities 
*/

class FlagPeriodicCounterparts_DataExchanger : public AOMD_DataExchanger
{
  int *n;
  list<mEntity *> *toSplit;
public :
  FlagPeriodicCounterparts_DataExchanger(int *nn, std::list<mEntity*> *l):n(nn),toSplit(l){}
  int tag() const override;
  void * AP_alloc_and_fill_buffer (mEntity *e, AOMD_SharedInfo &si, int) override;
  void receiveData (int pid, void *buf) override;   
};

int FlagPeriodicCounterparts_DataExchanger :: tag () const
{
  return 114;
}

void * FlagPeriodicCounterparts_DataExchanger :: AP_alloc_and_fill_buffer (mEntity *e, 
							     AOMD_SharedInfo &si, int d_tag)
{
  if(si.isPeriodic() && e->getData(AOMD_Util::Instance()->getAtt1()))
    {

      void *buf=malloc(sizeof(rp_int));

      rp_int *rp = (rp_int*)buf;
      rp->entity =  si.getRemotePointer();
      rp->i = 0; // dummy
      return buf;
    }
  else return nullptr;
}


void FlagPeriodicCounterparts_DataExchanger :: receiveData (int from, void *buf)
{
  rp_int *castbuf= (rp_int *)buf;
  if(!castbuf->entity->getData(AOMD_Util::Instance()->getAtt1()))
    {
      castbuf->entity->attachInt(AOMD_Util::Instance()->getAtt1(),++(*n));
      toSplit->push_back(castbuf->entity);
    }
}


/*-----------------------------------------------------------------*/
/*
  Main difference betweeen serial and parallel version of
  automatic refinement is that vertices id's in different
  processors must be different if they refer to different
  vertices and MUST be equal if they point to same vertices.
  This is not really important until we want to use load 
  balancing. Dynamic load balancing after adaptation will 
  become much more easy this way.

  We need 3 stages :

    -) first, we must count how many new vertices are needed 
    in each processors, this info must be broadcasted (first
    round of communication).
    -) second, we must insure that entities that are on partition
    boundaries and that are splitted will generate same new vertices 
    id's. As a rule of thumb, if both parts have decided to split, 
    then, the new vertex id will be the one of the smallest process
    rank(). This needs a second round of communications.
    -) Then, we split entities, using the id's provided.
*/


/*-------------------------------------------------------------------
This Callback is made to ensure that 
new vertices due to splitting have the 
same idea on different processors. 
If several processors have decided
to split the entity, the smallest
processor rank decides on vertex id.

  The strategy is as follow :
     -) An entity common to several processors may be splitted, 
        we attach on each processor an unique vertex id as a 
	predictor
     -) Then, the entity pass the information to each other processor
        where the entity is remotely defined and the id for the 
	new vertex will be the one from the smallest id which was generated
     -) This insures that vertex id's are identicla on different processors,
        this allows aomd to connect higher dimensional entities.
*/

int mMesh::split(list<mEntity*> &toSplit, mSplitCallbacks &f)
{
  // ----- 1 --------------------------------------------------------
  // globalMaxId is the highest vertex id in the whole distrib. mesh
  int localMaxId = theIdGenerator.getMaxValue();
  int startingVertex = localMaxId;


  //-----------------------------------------------------------------

  // ----- 2 --------------------------------------------------------
  // nbNewVertices is the number of new vertices in this partition

  //  ParUtil::Instance()->Msg(ParUtil::INFO,"s = %d\n",startingVertex);

  int nbNewVertices = 0;
  for(std::list<mEntity*>::const_iterator it = toSplit.begin();
      it != toSplit.end();
      ++it)
    {
      mEntity *e = *it;
      flagEntity(e,nbNewVertices,true,theOwnerManager);
    }

  //  ParUtil::Instance()->Msg(ParUtil::INFO,"entities flagged for splitting\n");
  //  printf("nbNewVertices = %d\n",nbNewVertices);
  FlagPeriodicCounterparts_DataExchanger fpc(&nbNewVertices,&toSplit);
  exchangeDataOnPartBdrys (fpc);
  //  printf("nbNewVertices = %d\n",nbNewVertices);

  //  ParUtil::Instance()->Msg(ParUtil::INFO,"periodic counterparts flagged\n");

  ParUtil::Instance()->Barrier(__LINE__,__FILE__);

  if(!toSplit.size())return 0;


  for(int i=1;i<4;i++)
    {
      for(iter it = begin(i) ; it != end(i) ; ++it)
	{
	  if((*it)->getData(AOMD_Util::Instance()->getAtt1()))
	    {	      
	      //             printf("re-flagging : ");(*it)->print();
	      setFlags(*it,startingVertex,true,theOwnerManager,theOwnerManager->isPeriodic(*it));
	    }
	}
    }

  //-----------------------------------------------------------------

  // ----- 4 --------------------------------------------------------
  // add to toSplit entities that are splitted on other processes
  
  //  ParUtil::Instance()->Msg(ParUtil::INFO,"counted %d new vertices\n",startingVertex);
  Refinement_DataExchanger de(&toSplit,false);
  exchangeDataOnPartBdrys (de);
  Refinement_DataExchanger dep(&toSplit,true);
  exchangeDataOnPartBdrys (dep);
  //  ParUtil::Instance()->Msg(ParUtil::INFO,"refinement data exchanger done\n");
  //-----------------------------------------------------------------

  // ----- 5 --------------------------------------------------------
  // do the split now
  splitgen(toSplit,f);
  //-----------------------------------------------------------------
  //  ParUtil::Instance()->Msg(ParUtil::INFO,"entities splitted\n");

  // ----- 6 --------------------------------------------------------
  // unflag entities
  {
  for(int i=1;i<4;i++)
    {
      for(iterall it = beginall(i) ; it != endall(i) ; ++it)
	{
	  if((*it)->getData(AOMD_Util::Instance()->getAtt1()))
	    {	      
	      setFlags(*it,startingVertex,false,theOwnerManager,theOwnerManager->isPeriodic(*it));
	    }
	}
    }
   }

  //  ParUtil::Instance()->Msg(ParUtil::INFO,"unflagging done\n");
  //-----------------------------------------------------------------

  // ----- 6 --------------------------------------------------------
  // recompute links between partitions
  bdryLinkSetup();

  //  printf("proc %d split done\n",ParUtil::Instance()->rank());

  return 1;
}

} // end of namespace
