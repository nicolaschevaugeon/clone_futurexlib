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
#include <fstream>
#include <cstdio>
#include <cstring>

#include "mAOMD.h"
#include "mMesh.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mFace.h"
#include "mTet.h"
#include "mHex.h"
#include "mPrism.h"

#include "mException.h"
#include "AOMD.h"

#ifndef SIM
#include "modeler.h"
#else
#include "SimModel.h"
#endif

using std::cout;
using std::ofstream;
using std::ostream;
using std::istream;
using std::endl;
using namespace std;
using namespace AOMD;
namespace AOMD {
/////
//   Function to clean Mesh from Gambit file. Merging the double node where
//   needed 
////
void cleandoubleVertex(mMesh *theMesh)
{

  printf("Cleaning and merging the double Vertex from Gambit File \n");
  printf(" This may take a while. A clean mesh file will be outputed:\n");
  printf(" ls.msh. Use it next time. \n");
  int kt=0;
  int kj=0;
  int todel,tokeep;
  std::map<int,int> todeltokeep;
  for(mMesh::iter it = theMesh->begin(0);it != theMesh->end(0); ++it)
  {
  kj++;
  mVertex *v1 = (mVertex*)(*it);
  int Id1=v1->getId();
  double v1x=v1->point()(0);
  double v1y=v1->point()(1);
  double v1z=v1->point()(2);

   for(mMesh::iter it2 = it;it2 != theMesh->end(0); ++it2)
   {
    mVertex *v2 = (mVertex*)(*it2);
    int Id2=v2->getId();
    double v2x=v2->point()(0);
    double v2y=v2->point()(1);
    double v2z=v2->point()(2);

    double d=(v1x-v2x)*(v1x-v2x)+(v1y-v2y)*(v1y-v2y)+(v1z-v2z)*(v1z-v2z);
    pGEntity gv1=v1->getClassification();
    pGEntity gv2=v2->getClassification();
    
    if (d<1.E-08&&Id1!=Id2) 
    {
     if (gv1==gv2||(GEN_type(gv1)>1&&GEN_type(gv2)>1&&gv1!=gv2))
     {
     kt++;
     todel=Id2; 
     tokeep=Id1;
     todeltokeep.insert(make_pair(todel,tokeep));
    }
    }
   }
  }
 
  theMesh->modifyState(0,1,true,0); // vertex -> edge
  theMesh->modifyState(0,2,true,0); // vertex -> face
  theMesh->modifyState(0,3,true,0); // vertex -> region

  std::map<int,int>::const_iterator itnode;
  std::set<mEntity *> todestroy;
  std::set<mEntity *> vertextodestroy;
  for(itnode=todeltokeep.begin();itnode!=todeltokeep.end();++itnode)
  {
   tokeep=(*itnode).second;
   todel=(*itnode).first;
   mEntity *v3 =theMesh->getVertex(todel);
   vertextodestroy.insert(v3);
   if (v3->isAdjacencyCreated(3))
   {
    for  (mAdjacencyContainer::iter irr= v3->begin(3);irr!= v3->end(3); ++irr)
   {
       mEntity *ent = (mEntity*)(*irr);
       todestroy.insert(ent);
   }
   }  
   if (v3->isAdjacencyCreated(2))
   {
   for  (mAdjacencyContainer::iter irr= v3->begin(2);irr!= v3->end(2); ++irr)
   {
       mEntity *ent = (mEntity*)(*irr);
       todestroy.insert(ent);
   }
   }
   
   if (v3->isAdjacencyCreated(1))
   {
    for  (mAdjacencyContainer::iter irr= v3->begin(1);irr!= v3->end(1); ++irr)
    {
       mEntity *ent = (mEntity*)(*irr);
       todestroy.insert(ent);
    }
   } 
}
  
  int k12=0;
  int s111=todestroy.size();
   printf("%d Entities to destroy \n",s111);
  std::set<mEntity*>::const_iterator itdest;
  for (itdest=todestroy.begin();itdest!=todestroy.end();++itdest)
  {
   k12++;  

   mEntity *ent2=(mEntity*)(*itdest);
       int gl=0;
       int Idtable[20];
       pGEntity gent=ent2->getClassification();
       for  (mAdjacencyContainer::iter in= ent2->begin(0);in!= ent2->end(0); ++in)
       {   
       //printf("1\n");

            mVertex *vert= (mVertex*)(*in);
            std::map<int, int>:: const_iterator ittodel = todeltokeep.find(vert->getId());
	    if (ittodel==todeltokeep.end()) 
	    {
	     Idtable[gl]=vert->getId();      
	     /// printf("20\n");

	    }
	    else{
	    ///printf("21\n");
	    Idtable[gl]=(*ittodel).second;}
	    gl++;
       }
       //EDGE,TRI,QUAD,HEX,PRISM,PYRAMID,TET,MIRROR,FACE
       if (ent2->getType()==mEntity::TET)	     
        {//printf("Tet %d %d %d %d\n", Idtable[0],Idtable[1],Idtable[2],Idtable[4]);

       theMesh->createTetWithVertices
        (theMesh->getVertex(Idtable[0]),
	 theMesh->getVertex(Idtable[1]), 
	 theMesh->getVertex(Idtable[2]),
	 theMesh->getVertex(Idtable[3]),
	 gent);}
       if (ent2->getType()==mEntity::TRI)
      theMesh->createFaceWithVertices
        (theMesh->getVertex(Idtable[0]),
	 theMesh->getVertex(Idtable[1]), 
	 theMesh->getVertex(Idtable[2]),
	 gent);
	if (ent2->getType()==mEntity::EDGE)
      theMesh->createEdge
        (theMesh->getVertex(Idtable[0]),
	 theMesh->getVertex(Idtable[1]),
	 gent); 
   
   
  //printf("%d deleted and reconstructed elements \n",k12);

   theMesh->DEL(ent2);
  }
  printf("%d deleted and reconstructed elements \n",k12);
  k12=0;
for (itdest=vertextodestroy.begin();itdest!=vertextodestroy.end();++itdest)
  {
   k12++;
   mEntity *ent2=(mEntity*)(*itdest); 
   theMesh->DEL(ent2);
  }
  
  printf("%d deleted nodes\n",k12);
  theMesh->modifyState(0,1,false,0); // vertex -> edge
  theMesh->modifyState(0,2,false,0); // vertex -> face
  theMesh->modifyState(0,3,false,0); // vertex -> face

  //  printf("Writing the cleaned mesh for futur use: cleaned.msh\n");

  //  AOMD::AOMD_Util::Instance()->ex_port("cleaned.msh",theMesh);
}


////////////////////////////




  void AOMD_Util::importGambitFile (const char *fName, mMesh *theMesh)
  {
   char line[256];
   int i;
   int NNOD, NELEM, NGRPS, NBSETS, NDFCD, NDFVL; 
   //NNOD=number of nodes
   //NELEM=number of elements
   //NGRPS=number of group
   //NBSETS=number of boundary set
   //NDFCD=dimension
   //NDFVL=
   int nodeID;
   double x,y,z;
   mVertex *vv;
   printf("import Gambit File \n");
   FILE *in = fopen (fName,"r");
   if(!in) printf("file %s not found",fName);
   if(!in)throw new mException (__LINE__,__FILE__,"unable to open file in loadGambitFile");
   
    do {
     fscanf(in,"%s",line);}
     while (strcmp(line,"NDFVL"));
     
     //read number of entity in the mesh
    fscanf(in,"%d %d %d %d %d %d",&NNOD, &NELEM, &NGRPS, &NBSETS, &NDFCD, &NDFVL );
    
    printf("reading nodes \n");   
     do {fscanf(in,"%s",line);} while (strcmp(line,"COORDINATES"));
     fscanf(in,"%s",line);
     if (NDFCD==3)
     {
     for (i=0;i<NNOD;i++)
      {
      fscanf(in,"%d %lf %lf %lf",&nodeID,&x,&y,&z);
      vv = theMesh->createVertex(nodeID,x,y,z,theMesh->getGEntity(1,3));
     };
     }
     else if (NDFCD==2)
     { 
      for(i=0;i<NNOD;i++)
      {
      fscanf(in,"%d %lf %lf",&nodeID,&x,&y);
      vv = theMesh->createVertex(nodeID,x,y,0.,theMesh->getGEntity(1,2));
      };
     }
     do {
     fscanf(in,"%s",line);
     }
     
     while (strcmp(line,"ELEMENTS/CELLS"));
     printf("reading elements \n");
     fscanf(in,"%s",line);
     std::map<int,mEntity *> e;
     int ivertex1,ivertex2,ivertex3,ivertex4;
     int ivertex5,ivertex6,ivertex7,ivertex8;
     int idelem,typelem,nnodes;
     for (i=0;i<NELEM;i++)
      {  
      fscanf(in,"%d %d %d",&idelem,&typelem,&nnodes);
      switch (typelem)
       {
       case 1 : //Edge
        {
	 fscanf(in,"%d %d",&ivertex1,&ivertex2);
         e[idelem]=theMesh->createEdge(theMesh->getVertex(ivertex1),theMesh->getVertex(ivertex2),theMesh->getGEntity(1,NDFCD));
	 break;
        }
       case 2 : //Quadrilateral
        {	
	 fscanf(in,"%d %d %d %d",&ivertex1,&ivertex2,&ivertex3,&ivertex4);
          e[idelem]=theMesh->createFaceWithVertices(theMesh->getVertex(ivertex1),theMesh->getVertex(ivertex2),
	 theMesh->getVertex(ivertex3),theMesh->getVertex(ivertex4),theMesh->getGEntity(1,NDFCD));
	 break;

        }
       case 3 : //Triangle
        {	 
	 fscanf(in,"%d %d %d",&ivertex1,&ivertex2,&ivertex3);
          e[idelem]=theMesh->createFaceWithVertices(theMesh->getVertex(ivertex1),theMesh->getVertex(ivertex2),
	 theMesh->getVertex(ivertex3),theMesh->getGEntity(1,NDFCD));
        break;
       }
       case 4 : //Brick (Hex)
        {
	 fscanf(in,"%d %d %d %d %d %d %d %d",&ivertex1,&ivertex2,&ivertex3,&ivertex4,
	  &ivertex5,&ivertex6,&ivertex7,&ivertex8);
         e[idelem]=theMesh->createHexWithVertices(theMesh->getVertex(ivertex1),theMesh->getVertex(ivertex2),
	  theMesh->getVertex(ivertex3),theMesh->getVertex(ivertex4),theMesh->getVertex(ivertex5),
	  theMesh->getVertex(ivertex6),theMesh->getVertex(ivertex7),theMesh->getVertex(ivertex8)
	  ,theMesh->getGEntity(1,NDFCD));
         break;
	}
       case 5 : //Wedge (Prism)
        {	 
	 fscanf(in,"%d %d %d %d %d %d",&ivertex1,&ivertex2,&ivertex3,&ivertex4,
	  &ivertex5,&ivertex6);
         e[idelem]=theMesh->createPrismWithVertices(theMesh->getVertex(ivertex1),theMesh->getVertex(ivertex2),
	  theMesh->getVertex(ivertex3),theMesh->getVertex(ivertex4),theMesh->getVertex(ivertex5),
	 theMesh->getVertex(ivertex6)
	  ,theMesh->getGEntity(1,NDFCD));
         break;
	}
       case 6 : //Tetrahedron
       {
	 fscanf(in,"%d %d %d %d",&ivertex1,&ivertex2,&ivertex3,&ivertex4);
         e[idelem]=theMesh->createTetWithVertices(theMesh->getVertex(ivertex1),theMesh->getVertex(ivertex2),
	  theMesh->getVertex(ivertex3),theMesh->getVertex(ivertex4)
	  ,theMesh->getGEntity(1,NDFCD));
         break;
	}
       case 7 : //Pyramid
        {         
	fscanf(in,"%d %d %d %d %d ",&ivertex1,&ivertex2,&ivertex3,&ivertex4,
	  &ivertex5);
	  printf("AOMD did not contains Pyramids element at the time mGambit.cc was programmed");
       // e[idelem]=(mEntity *)theMesh->createPyramidWithVertices(theMesh->getVertex(ivertex1),theMesh->getVertex(ivertex2),
	//  theMesh->getVertex(ivertex3),theMesh->getVertex(ivertex4),theMesh->getVertex(ivertex5),
	 // theMesh->getGEntity(0,0));
	 break;
	}
     }
     }
   
     
   
   if (NDFCD==2) theMesh->modifyState(0,2,true,0); // vertex -> face
   if (NDFCD==3) theMesh->modifyState(0,3,true,0); // vertex -> region
   int NELEMGROUP,dump1,NFLAGS,idGROUP,k,elemId;
   int FLAG[10];
    int verif=0;
    printf("reading Groups \n");
    for (i=0;i<NGRPS;i++)
   {
    do {fscanf(in,"%s",line);}
     while (strcmp(line,"GROUP"));
     fscanf(in,"%s",line);
//     fscanf(in,"GROUP: %d ELEMENTS: %d MATERIAL: %d NLFAGS: %10d",&idGROUP,&NELEMGROUP,&dump1,&dump2);
          fscanf(in," GROUP: %d ELEMENTS: %d MATERIAL: %d NFLAGS: %d",&idGROUP,&NELEMGROUP,&dump1,&NFLAGS);
     fscanf(in,"%s",line);
     for (k=0;k<NFLAGS;k++) fscanf(in,"%d",&FLAG[k]);
     for (k=0;k<NELEMGROUP;k++)
      { 
        verif=fscanf(in,"%d",&elemId);
        e[elemId]->classify(theMesh->getGEntity(idGROUP,NDFCD));
	for  (mAdjacencyContainer::iter irr= e[elemId]->begin(0);irr!= e[elemId]->end(0); ++irr)
	{
	   mEntity *ent = (mEntity*)(*irr);
	   ent->classify(theMesh->getGEntity(idGROUP,NDFCD));
	}
	
      }
   }

   int count=0;
   printf("reading BoundaryConditions \n");
   char NAME[256];
   int ITYPE, NENTRY, NVALUES, IBCODE1;
   int ELEMID, ELEMTYPE, FACE;
   mEntity *FACEe;
   mVertex *vertex0,*vertex1,*vertex2,*vertex3;
    do
    {do {verif=fscanf(in,"%s",line);}
       while (strcmp(line,"CONDITIONS") && verif !=EOF);
       if (verif != EOF) 
       {
       count++;
       fscanf(in,"%s",line);
       fscanf(in,"%s %d %d %d %d",NAME, &ITYPE, &NENTRY, &NVALUES, &IBCODE1);
       if (ITYPE==1)
        {
        for (i=0;i<NENTRY;i++)
         {
          fscanf(in,"%d %d %d",&ELEMID, &ELEMTYPE, &FACE);    
          FACEe=e[ELEMID]->getTemplate(FACE-1,NDFCD-1,0);
          // !!!!!!!! GAMBIT NUMBER FACES FROM 1 to 4 WHILE AOMD NUMBER FROM 0 to 3 (FACE-1) 
          switch (ELEMTYPE)
           {
            case 2 : //Quadrilateral
	     {
	      vertex0=(mVertex *)FACEe->get(0,0);
              vertex1=(mVertex *)FACEe->get(0,1);
              FACEe=theMesh->createEdge(vertex0,vertex1,theMesh->getGEntity(count*10,1));
	     for  (mAdjacencyContainer::iter irr= FACEe->begin(0);irr!= FACEe->end(0); ++irr)
	     {
	       mEntity *ent = (mEntity*)(*irr);
	       ent->classify(theMesh->getGEntity(count*10,1));
	     }
	
	     
	     break;
	    }
	   case 3 : //Triangle
	    {
	     vertex0=(mVertex *)FACEe->get(0,0);
             vertex1=(mVertex *)FACEe->get(0,1);
             FACEe=theMesh->createEdge(vertex0,vertex1,theMesh->getGEntity(count*10,1));
	     for  (mAdjacencyContainer::iter irr= FACEe->begin(0);irr!= FACEe->end(0); ++irr)
	     {
	       mEntity *ent = (mEntity*)(*irr);
	       ent->classify(theMesh->getGEntity(count*10,1));
	     }
	     break;
	    }
       	   case 4 : //Brick (Hex)
	    {
	     vertex0=(mVertex *)FACEe->get(0,0);
             vertex1=(mVertex *)FACEe->get(0,1);
             vertex2=(mVertex *)FACEe->get(0,2);
	     vertex3=(mVertex *)FACEe->get(0,3);
	     FACEe=theMesh->createFaceWithVertices(vertex0,vertex1,vertex2,vertex3,theMesh->getGEntity(count*10,2));
	     break;
	    }
           case 5 : //Wedge (Prism)
            {
	     switch (FACE-1)
	      {
	       case 0 :
	       case 1 :
	       case 2 :
	        vertex0=(mVertex *)FACEe->get(0,0);
                vertex1=(mVertex *)FACEe->get(0,1);
                vertex2=(mVertex *)FACEe->get(0,2);
	        vertex3=(mVertex*)FACEe->get(0,3);
	        FACEe=theMesh->createFaceWithVertices(vertex0,vertex1,vertex2,vertex3,theMesh->getGEntity(count*10,2));
	        break;
	       case 3:
	       case 4:
	        vertex0=(mVertex *)FACEe->get(0,0);
                vertex1=(mVertex *)FACEe->get(0,1);
                vertex2=(mVertex *)FACEe->get(0,2);	      
	        FACEe=theMesh->createFaceWithVertices(vertex0,vertex1,vertex2,theMesh->getGEntity(count*10,2));
	        break;
	      }
	      break;
	    }
	   case 6 : //Tetrahedron
	    {
	     vertex0=(mVertex *)FACEe->get(0,0);
             vertex1=(mVertex *)FACEe->get(0,1);
             vertex2=(mVertex *)FACEe->get(0,2);
             FACEe=theMesh->createFaceWithVertices(vertex0,vertex1,vertex2,theMesh->getGEntity(count*10,2));
	     break;
	    }
           case 7 : //Pyramid
            { 
	     printf("AOMD did not contains Pyramids element at the time mGambit.cc was programmed");
             break;
	    }
       }
      }
      }//
      }
     }
     while (verif != EOF);
   if(in) fclose(in);

   printf("end import Gambit File\n");
     
     //AOMD::AOMD_Util::Instance()->exportGmshFile("avcleaned.msh",theMesh);

   // cleandoubleVertex(theMesh);
  }






}



