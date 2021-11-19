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
#include "mAOMD.h"
#include "AOMD_OwnerManager.h"
#include "ParUtil.h"
#include "mEdge.h"
#include "mException.h"
#include "mFace.h"
#include "mMesh.h"
#include "mMirrorEntity.h"
#include "mTet.h"
#include "mVertex.h"

#include "modeler.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>

using std::cout;
using std::endl;
using std::istream;
using std::ofstream;
using std::ostream;
using std::vector;
using namespace std;
namespace AOMD
{
// AOMD_Util  AOMD_Util::instance;
AOMD_Util::object_creator AOMD_Util::create_object;

AOMD_Util::AOMD_Util()
{
   _ATTPARAMETRIC = lookupMeshDataId("_parametric");
   _ATTPARENT = lookupMeshDataId("_parent");
   _ATTFMOD = lookupMeshDataId("_fmod");
   _ATTEMOD = lookupMeshDataId("_emod");
   _ATTID = lookupMeshDataId("_id");
   _ATTDN = lookupMeshDataId("_dn");
   _ATTSIZE = lookupMeshDataId("_size");
   _ATT1 = lookupMeshDataId("_ATT1");
   _ATT2 = lookupMeshDataId("_ATT2");
   _ATT3 = lookupMeshDataId("_ATT3");
   _ATTWEIGHT = lookupMeshDataId("_weight");
   currentProcessor = -1;
}

AOMD_Util::~AOMD_Util()
{
   //    std::cout << "DETROYING AOMD UTIL" << std::endl;
}

AOMD_Util *AOMD_Util::Instance()
{
   static AOMD_Util instance;
   create_object.do_nothing();
   return &instance;
}

void classifyUnclassifiedVerices(mMesh *m)
{
   for (int i = 1; i < 4; i++)
   {
      // int k  = 0;
      if (m->size(i))
      {
         for (mMesh::iterall it = m->beginall(i); it != m->endall(i); ++it)
         {
            mEntity *e = *it;
            //	  printf("----%d(%p)---\n",k++,e);
            //	  e->print();
            for (int j = 0; j < e->size(0); j++)
            {
               mEntity *v = e->get(0, j);
               //	      v->print();
               if (!v->getClassification() && e->getClassification()) v->classify(e->getClassification());
               if (!e->getClassification()) e->print();
            }
         }
      }
   }
}

void AOMD_Util::import(const char *fName, mMesh *theMesh)
{
   char ext[6];
   strcpy(ext, fName + (strlen(fName) - 4));
   if (!strcmp(ext, ".sms"))
   {
      importSmsFile(fName, theMesh);
   }
   else if (!strcmp(ext, ".msh"))
   {
      importDGFile(fName, theMesh);
   }
   else if (!strcmp(ext, ".pat"))
   {
      char text[256];
      sprintf(text, "\".pat\" import format not done yet, cannot import %s", fName);
      throw mException(__LINE__, __FILE__, text);
   }
   else if (!strcmp(ext, ".neu"))
   {
      importGambitFile(fName, theMesh);
   }
   else
   {
      char text[256];
      sprintf(text, "unknown extension %s in file %s", ext, fName);
      throw new mException(__LINE__, __FILE__, text);
   }
   BuildGTopology(theMesh);
   //  print_topology ("topo.dat",theMesh);
}

void AOMD_Util::import_oneLevel(const char *fName, mMesh *theMesh)
{
   char ext[6];
   strcpy(ext, fName + (strlen(fName) - 4));
   if (!strcmp(ext, ".sms"))
   {
#ifndef SIM
      importSmsFile_oneLevel(fName, theMesh);
#endif
   }
   else if (!strcmp(ext, ".msh"))
   {
      //	importDGFile_oneLevel(fName, theMesh);
      char text[256];
      sprintf(text, "\".pat\" import format not done yet, cannot import %s", fName);
      throw mException(__LINE__, __FILE__, text);
   }
   else if (!strcmp(ext, ".pat"))
   {
      char text[256];
      sprintf(text, "\".pat\" import format not done yet, cannot import %s", fName);
      throw mException(__LINE__, __FILE__, text);
   }
   else if (!strcmp(ext, ".neu"))
   {
      char text[256];
      sprintf(text, "\".pat\" import format not done yet, cannot import %s", fName);
      throw mException(__LINE__, __FILE__, text);
   }
   else
   {
      char text[256];
      sprintf(text, "unknown extension %s in file %s", ext, fName);
      throw new mException(__LINE__, __FILE__, text);
   }
   //   BuildGTopology(theMesh);
   //  print_topology ("topo.dat",theMesh);
}

extern bool isEntityOnCurrentProcessor(mEntity *e, int currentProcessor);

void AOMD_Util::ex_port(const char *fName, const mMesh *theMesh, bool reduce_to_minimum)
{
   char ext[6];
   char text[256];
   strcpy(ext, fName + (strlen(fName) - 4));
   if (!strcmp(ext, ".sms"))
   {
      // sprintf(text,"\".sms\" export format not done yet, cannot export %s",fName);
      // throw mException (__LINE__,__FILE__,text);
      exportSmsFile(fName, theMesh);
   }
   else if (!strcmp(ext, ".msh"))
   {
      ofstream ofs(fName);
      exportDGFile(ofs, theMesh, reduce_to_minimum);
      ofs.close();
   }
   else if (!strcmp(ext, ".pat"))
   {
      char text[256];
      sprintf(text, "\".pat\" export format not done yet, cannot export %s", fName);
      throw mException(__LINE__, __FILE__, text);
      //      exportPatranNeutralFile(fName, theMesh);
   }
   else if (!strcmp(ext, ".geo"))
   {
      exportEnsightFile(fName, theMesh, reduce_to_minimum);
      //      export Ensight geometry file
      // exportEnsightFile(fName, theMesh);
   }
   else if (!strcmp(ext, ".vtu"))
   {
      exportVTUFile(fName, theMesh, reduce_to_minimum);
   }
   else
   {
      sprintf(text, "unknown extension %s in file %s", ext, fName);
      throw mException(__LINE__, __FILE__, text);
   }
}

void AOMD_Util::importDGFile(const char *fName, mMesh *theMesh)
{
   std::ifstream ifs(fName);
   importDGFile(ifs, theMesh);
   ifs.close();
}

void printGmshFace(ostream &fs, mFace *f, int &k, int recur)
{
   //  if(recur == 2)return;

   mVertex *vx[4];
   for (int i = 0; i < f->size(0); i++)
   {
      mVertex *v = (mVertex *)f->get(0, i);
      vx[i] = v;
   }

   {
      for (int i = 0; i < recur; i++) fs << " ";
   }
   int nbrec = (f->isAdjacencyCreated(2)) ? f->size(2) : 0;
   int typ = (f->size(0) == 3) ? 2 : 3;
   if (typ == 2)
      fs << k++ << " " << typ << " " << GEN_tag(f->getClassification()) << " " << GEN_tag(f->getClassification()) << " " <<
          // -nbrec << " " <<
          3 << " " << vx[0]->getId() << " " << vx[1]->getId() << " " << vx[2]->getId() << "\n";

   else
      fs << k++ << " " << typ << " " << GEN_tag(f->getClassification()) << " " << -nbrec << " " << 4 << " " << vx[0]->getId()
         << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << "\n";

   {
      for (int i = 0; i < nbrec; i++) printGmshFace(fs, (mFace *)f->get(2, i), k, recur + 1);
   }
}

void printGmshRegion(ostream &fs, mRegion *f, int &k, int recur)
{
   //  if(recur == 2)return;

   mVertex *vx[8];

   for (int i = 0; i < f->size(0); i++)
   {
      mVertex *v = (mVertex *)f->get(0, i);
      vx[i] = v;
   }

   {
      for (int i = 0; i < recur; i++) fs << " ";
   }
   int nbrec = (f->isAdjacencyCreated(3)) ? f->size(3) : 0;
   switch (f->getType())
   {
      case mEntity::PRISM:
         fs << k++ << " " << 6 << " " << GEN_tag(f->getClassification()) << " " << -nbrec << " " << 6 << " " << vx[0]->getId()
            << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << " " << vx[4]->getId() << " "
            << vx[5]->getId() << "\n";
         break;
      case mEntity::HEX:
         fs << k++ << " " << 5 << " " << GEN_tag(f->getClassification()) << " " << -nbrec << " " << 8 << " " << vx[0]->getId()
            << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << " " << vx[4]->getId() << " "
            << vx[5]->getId() << " " << vx[6]->getId() << " " << vx[7]->getId() << "\n";
         break;
      case mEntity::TET:
         fs << k++ << " " << 4 << " " << GEN_tag(f->getClassification()) << " " << GEN_tag(f->getClassification()) << " " <<
             //	  -nbrec << " " <<
             4 << " " << vx[0]->getId() << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << "\n";
         break;
      default:
      {
         std::cout << "Error in file " << __FILE__ << ":" << __LINE__ << " entity of type " << f->getType()
                   << " not handled in switch." << std::endl;
         throw;
         break;
      }
   }
   for (int i = 0; i < nbrec; i++) printGmshRegion(fs, (mRegion *)f->get(3, i), k, recur + 1);
}

void printGmshEdge(ostream &fs, mEdge *f, int &k, int recur)
{
   //  if(recur == 2)return;
   mVertex *vx[2];
   for (int i = 0; i < f->size(0); i++)
   {
      mVertex *v = (mVertex *)f->get(0, i);
      vx[i] = v;
   }
   for (int i = 0; i < recur; i++) fs << " ";
   int nbrec = (f->isAdjacencyCreated(1)) ? f->size(1) : 0;
   fs << k++ << " 1 " << GEN_tag(f->getClassification()) << " " << -nbrec << " " << 2 << " " << vx[0]->getId() << " "
      << vx[1]->getId() << "\n";

   for (int i = 0; i < nbrec; i++) printGmshEdge(fs, (mEdge *)f->get(1, i), k, recur + 1);
}

void AOMD_Util::importMichelGrid(istream &input, mMesh *mesh)
{
   int nVerts, nTets, nTris;
   input >> nVerts >> nTets >> nTris;

   cout << "o Grid sizes: " << nVerts << " " << nTets << " " << nTris << endl;

   cout << "o Reading vertices ... " << endl;

   int iVert;
   double x, y, z;
   pGEntity geoEntity = nullptr;
   for (iVert = 0; iVert < nVerts; iVert++)
   {
      input >> x >> y >> z;
      mesh->createVertex(iVert, x, y, z, geoEntity);
   }

   cout << "o Reading Tets ... " << endl;

   int iTet;
   int iD0, iD1, iD2, iD3;
   for (iTet = 0; iTet < nTets; iTet++)
   {
      input >> iD0 >> iD1 >> iD2 >> iD3;
      iD0--;
      iD1--;
      iD2--;
      iD3--;
      mesh->createTetWithVertices(iD0, iD1, iD2, iD3, geoEntity);
   }
   cout << "o Reading done ... " << endl;
}

void AOMD_Util::exportGmshFile(const char *fName, const mMesh *theMesh)
{
   ofstream f(fName);
   exportGmshFile(f, theMesh);
   f.close();
}

void AOMD_Util::exportGmshFile(ostream &f, const mMesh *theMesh)
{
   f << "$NOE\n";
   f << theMesh->size(0) << "\n";
   for (mMesh::iter it = theMesh->begin(0); it != theMesh->end(0); ++it)
   {
      mVertex *v = (mVertex *)(*it);
      f << v->getId() << " " << v->point()(0) << " " << v->point()(1) << " " << v->point()(2) << "\n";
   }
   f << "$ENDNOE\n";
   int k = 0;

   int lev = 3;
   if (!theMesh->size(3)) lev = 2;
   if (!theMesh->size(2) && !theMesh->size(3)) lev = 1;

   //  printf("exporting mesh %d with lev %d %d\n",theMesh->getId(),lev,theMesh->size(1));
   /*  {
       for(mMesh::iter it = theMesh->begin(0);it != theMesh->end(0); ++it)
       {
       mEdge *edge = (mEdge*)(*it);
       if(lev==1 || edge->getClassification()->getLevel() == 1)
       k++;
       }
       }
   */
   {
      for (mMesh::iter it = theMesh->begin(1); it != theMesh->end(1); ++it)
      {
         mEdge *edge = (mEdge *)(*it);
         if (lev == 1 || GEN_type(edge->getClassification()) == 1 || GEN_tag(edge->getClassification()) < 0) k++;
      }
   }
   {
      for (mMesh::iterall it = theMesh->beginall(2); it != theMesh->endall(2); ++it)
      {
         mFace *face = (mFace *)(*it);
         if (lev == 2 || GEN_type(face->getClassification()) == 2 || GEN_tag(face->getClassification()) < 0) k++;
      }
   }

   f << "$ELM\n";
   f << k + theMesh->size(3) << "\n";

   k = 1;
   {
      for (mMesh::iter it = theMesh->begin(1); it != theMesh->end(1); ++it)
      {
         mEdge *edge = (mEdge *)(*it);
         if (lev == 1 || GEN_type(edge->getClassification()) == 1) printGmshEdge(f, edge, k, 0);
      }
   }
   {
      for (mMesh::iterall it = theMesh->beginall(2); it != theMesh->endall(2); ++it)
      {
         mFace *face = (mFace *)(*it);
         if (lev == 2 || GEN_type(face->getClassification()) == 2) printGmshFace(f, face, k, 0);
      }
   }
   {
      for (mMesh::iter it = theMesh->begin(3); it != theMesh->end(3); ++it)
      {
         printGmshRegion(f, (mRegion *)(*it), k, 0);
      }
   }
   f << "$ENDELM\n";
   f << "$FIN\n";
}

void AOMD_Util::exportDGFile(ostream &f, const mMesh *theMesh, bool reduce_to_minimum)
{
   //-- N O D   S E C T I O N  --------------------------------------------------
   // exporting vertices
   // format : id x y z
   f << "$NOE\n";
   int NPARAM = 0;
   int NV = 0;
   int NV_GV = 0;
   int NV_MIRROR = 0;
   for (mMesh::iterall it = theMesh->beginall(0); it != theMesh->endall(0); ++it)
   {
      if (isEntityOnCurrentProcessor((*it), currentProcessor))
      {
         mVertex *v = (mVertex *)(*it);
         if (v->getType() == mEntity::VERTEX)
            NV++;
         else if (v->getType() == mEntity::MIRROR)
            NV_MIRROR++;
         // buggy else if is now if
         if (GEN_type(v->getClassification()) == 0 || v->getCommonBdry()) NV_GV++;
         if (v->getData(getParametric())) NPARAM++;
      }
   }

   f << NV << "\n";
   {
      for (mMesh::iterall it = theMesh->beginall(0); it != theMesh->endall(0); ++it)
      {
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            mVertex *v = (mVertex *)(*it);
            if (v->getType() == mEntity::VERTEX)
               f << v->getId() << " " << v->point()(0) << " " << v->point()(1) << " " << v->point()(2) << "\n";
         }
      }
   }
   f << "$ENDNOE\n";

   //-- S Y M   S E C T I O N --------------------------------------------------
   // exporting vertex associations i.e. mirror vertices : list of equivalent
   // vertices

   if (NV != theMesh->size(0))
   {
      f << "$SYM\n";
      f << NV_MIRROR << "\n";
      for (mMesh::iterall it = theMesh->beginall(0); it != theMesh->endall(0); ++it)
      {
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            mVertex *v = (mVertex *)(*it);
            if (v->getType() == mEntity::MIRROR)
            {
               mMirrorVertex *v = (mMirrorVertex *)(*it);
               f << v->nbCopies() << " ";
               for (int i = 0; i < v->nbCopies(); i++) f << v->getCopy(i)->getId() << " ";
               f << "\n";
            }
         }
      }
      f << "$ENDSYM\n";
   }

   //-- P R L   S E C T I O N --------------------------------------------------
   // exporting parallel associations

   f << "$PRL\n";
   int N(0);
   for (AOMD_OwnerManager::iter_top itx = theMesh->theOwnerManager->begin_top(); itx != theMesh->theOwnerManager->end_top();
        ++itx)
      N++;
   f << N << "\n";
   {
      for (AOMD_OwnerManager::iter_top itx = theMesh->theOwnerManager->begin_top(); itx != theMesh->theOwnerManager->end_top();
           ++itx)
      {
         int id = (*itx).first;
         int to = (*itx).second;
         f << id << " " << to << std::endl;
      }
   }
   f << "$ENDPRL\n";

   //-- P R M   S E C T I O N --------------------------------------------------
   // exporting parametric coordinates of vertices

   if (NPARAM)
   {
      f << "$PRM\n";
      f << NPARAM << "\n";
      for (mMesh::iterall it = theMesh->beginall(0); it != theMesh->endall(0); ++it)
      {
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            mVertex *v = (mVertex *)(*it);
            if (v->getData(getParametric()))
            {
               Trellis_Util::mVector par = v->getAttachedVector(getParametric());
               f << v->getId() << " " << par(0) << " " << par(1) << " " << par(2) << "\n";
            }
         }
      }
      f << "$ENDPRM\n";
   }

   //-- E L M   S E C T I O N --------------------------------------------------
   // exporting elements

   int NBELM = NV_GV;
   for (int dim = 1; dim <= 3; dim++)
      for (mMesh::iterall it = theMesh->beginall(dim); it != theMesh->endall(dim); ++it)
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            int pID_dim = 4;
            if ((*it)->getCommonBdry()) pID_dim = (*it)->getCommonBdry()->getLevel();

            if (!reduce_to_minimum || GEN_type((*it)->getClassification()) < dim + 1 || pID_dim < dim + 1) NBELM++;
         }
   f << "$ELM\n";
   f << NBELM << "\n";
   int k = 1;
   {
      for (mMesh::iterall it = theMesh->beginall(0); it != theMesh->endall(0); ++it)
      {
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            /// Vertex classified on a model vertex or on
            /// a parallel model vertex
            mVertex *vertex = (mVertex *)(*it);
            int pID = 0;
            if (vertex->getCommonBdry()) pID = vertex->getCommonBdry()->getId();

            if (GEN_type(vertex->getClassification()) < 1 || vertex->getCommonBdry())
               f << k++ << " " << 15 << " " << GEN_tag(vertex->getClassification()) << " " << pID << " 1 " << vertex->getId()
                 << "\n";
         }
      }
   }
   {
      for (mMesh::iterall it = theMesh->beginall(1); it != theMesh->endall(1); ++it)
      {
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            mEdge *edge = (mEdge *)(*it);

            int pID_dim = 4;
            if (edge->getCommonBdry()) pID_dim = edge->getCommonBdry()->getLevel();

            if (!reduce_to_minimum || GEN_type(edge->getClassification()) < 2 || pID_dim < 2)
               if (!edge->parent()) printGmshEdge(f, edge, k, 0);
         }
      }
   }
   {
      for (mMesh::iterall it = theMesh->beginall(2); it != theMesh->endall(2); ++it)
      {
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            mFace *face = (mFace *)(*it);

            int pID_dim = 4;
            if (face->getCommonBdry()) pID_dim = face->getCommonBdry()->getLevel();

            if (!reduce_to_minimum || GEN_type(face->getClassification()) < 3 || pID_dim < 3)
               if (!face->parent()) printGmshFace(f, face, k, 0);
         }
      }
   }
   {
      for (mMesh::iterall it = theMesh->beginall(3); it != theMesh->endall(3); ++it)
      {
         if (isEntityOnCurrentProcessor((*it), currentProcessor))
         {
            mRegion *region = (mRegion *)(*it);
            if (!region->parent()) printGmshRegion(f, region, k, 0);
         }
      }
   }
   f << "$ENDELM\n";
}

void AOMD_Util::exportEnsightFile(const char *fileName, const mMesh *theMesh, bool reduce_to_minimum)
{
   printf("AOMD Warning: exportEnsightFile is not implemented yet\n");
   // exit;

   /*
       FILE *f = fopen (fileName,"w");


       //-- N O D   S E C T I O N  --------------------------------------------------
       // exporting vertices
       // format : id x y z

       int NV = 0;
       int NV_GV = 0;
       int maxdim=theMesh->getDim();
       int countfluide=1;
       int countboundary=1;
       int countsupport=1;
       for(mMesh::iter it = theMesh->begin(0);it != theMesh->end(0); ++it)
         {
           mVertex *v = (mVertex*)(*it);
           if(v->getType() == mEntity::VERTEX)NV++;
           //buggy else if is now if
           if(v->getClassification()->type() == 0)
             NV_GV++;
         }
     fprintf(f,"Problem geometrie\n\nnode id given\nelement id assign\ncoordinates\n%8d\n", NV);

       { for(mMesh::iter it = theMesh->begin(0);it != theMesh->end(0); ++it)
         {
           mVertex *v = (mVertex*)(*it);
           if(v->getType() == mEntity::VERTEX)
             fprintf(f,"%8d%12.5e%12.5e%12.5e\n",(v->getId()+1),v->point()(0),v->point()(1),v->point()(2)); //add one to node ID
   because Ensight don't support node Id less than one.
         }
       }

       //-- E L M   S E C T I O N --------------------------------------------------
       // exporting elements
       int NBELM = NV_GV;
       int NBVERTEX=0;
       int NBEDGE=0;
       int NBTRI=0;
       int NBQUAD=0;
       int NBPYRAMID=0;
       int NBHEX=0;
       int NBPRISM=0;
       int NBTET=0;

   //    for (mMesh::giter git=theMesh->begin(); git != theMesh->end();++git)

       {
       AOMD::gEntity *thegEntity=*(git);
       int NBELM = NV_GV;
       int NBVERTEX=0;
       int NBEDGE=0;
       int NBTRI=0;
       int NBQUAD=0;
       int NBPYRAMID=0;
       int NBHEX=0;
       int NBPRISM=0;
       int NBTET=0;
       int EntityId=thegEntity->getId();
       int dim=thegEntity->getLevel();

       for(mMesh::iter it = theMesh->begin(dim);it != theMesh->end(dim); ++it)
       {
           AOMD::mEntity *meshEntity=*(it);
           AOMD::gEntity *geoEntity= meshEntity->getClassification();
           if (geoEntity==thegEntity)
           //printf("yes");
           {
            mEntity::mType TYPELEM=meshEntity->getType ();
            if(TYPELEM == mEntity::VERTEX) NBVERTEX++;
            if(TYPELEM == mEntity::EDGE) NBEDGE++;
            if(TYPELEM == mEntity::TRI) NBTRI++;
            if(TYPELEM == mEntity::QUAD) NBQUAD++;
            if(TYPELEM == mEntity::TET) NBTET++;
            if(TYPELEM == mEntity::HEX) NBHEX++;
            if(TYPELEM == mEntity::PRISM) NBPRISM++;
           }
       }
       if (dim==maxdim)
        {
          fprintf(f,"part %d%d \nFluide_domaine%d\n",dim,EntityId,countfluide);
          countfluide++;
        }
       else if (dim==maxdim-1)
        {
         fprintf(f,"part %d%d \nBoundary%d\n",dim,EntityId,countboundary);
         countboundary++;
        }
       else
        {
          fprintf(f,"part %d%d \nGeometric_support%d\n",dim,EntityId,countsupport);
          countsupport++;
        }
       if (NBVERTEX!=0)
       {
          fprintf(f,"point \n%8d\n",NBVERTEX);
          for(mMesh::iter it = theMesh->begin(0);it != theMesh->end(0); ++it)
          {
             AOMD::mEntity *meshEntity=*(it);
             AOMD::gEntity *geoEntity= meshEntity->getClassification();
             if (geoEntity==thegEntity)
              {
               mEntity::mType TYPELEM=meshEntity->getType ();
               if(TYPELEM == mEntity::VERTEX)
               {
                   fprintf(f,"%8d\n",meshEntity->getId()+1);
               }
              }
           }
         }
       if (NBEDGE!=0)
       {
       fprintf(f,"bar2 \n%8d\n",NBEDGE);
       for(mMesh::iter it = theMesh->begin(1);it != theMesh->end(1); ++it)
         {
           AOMD::mEntity *meshEntity=*(it);
           mEntity::mType TYPELEM=meshEntity->getType ();
           AOMD::gEntity *geoEntity= meshEntity->getClassification();
           if (geoEntity==thegEntity)
             {
             if(TYPELEM == mEntity::EDGE)
               {
                   for (int i=0;i<2;i++)
                   {
                           mVertex *V=(mVertex *)meshEntity->get(0,i);
                           fprintf(f,"%8d",(V->getId()+1));
                   }
                   fprintf(f,"\n");
                }
               }
           }
        }

       if (NBTRI!=0)
       {
       fprintf(f,"tria3 \n%8d\n",NBTRI);
       for(mMesh::iter it = theMesh->begin(2);it != theMesh->end(2); ++it)
        {
           AOMD::mEntity *meshEntity=*(it);
           AOMD::gEntity *geoEntity= meshEntity->getClassification();
           if (geoEntity==thegEntity)
             {
             mEntity::mType TYPELEM=meshEntity->getType ();
             if(TYPELEM == mEntity::TRI)
               {
                   for (int i=0;i<3;i++)
                   {
                           mVertex *V=(mVertex *)meshEntity->get(0,i);
                           fprintf(f,"%8d",(V->getId()+1));
                   }
                   fprintf(f,"\n");
               }
             }
        }
       }
       if (NBQUAD!=0)
       {
       fprintf(f,"quad4 \n%8d\n",NBQUAD);
       for(mMesh::iter it = theMesh->begin(2);it != theMesh->end(2); ++it)
         {
           AOMD::mEntity *meshEntity=*(it);
           AOMD::gEntity *geoEntity= meshEntity->getClassification();
             if (geoEntity==thegEntity)
             {
              mEntity::mType TYPELEM=meshEntity->getType ();
              if(TYPELEM == mEntity::QUAD)
               {
               for (int i=0;i<4;i++)
                   {
                           mVertex *V=(mVertex *)meshEntity->get(0,i);
                           fprintf(f,"%8d",(V->getId()+1));
                   }
                   fprintf(f,"\n");
                }
             }
          }
        }

       if (NBTET!=0)
       {
       fprintf(f,"tetra4 \n%8d\n",NBTET);
       for(mMesh::iter it = theMesh->begin(3);it != theMesh->end(3); ++it)
         {
           AOMD::mEntity *meshEntity=*(it);
           mEntity::mType TYPELEM=meshEntity->getType ();
           AOMD::gEntity *geoEntity= meshEntity->getClassification();
           if (geoEntity==thegEntity)
           {
            if(TYPELEM == mEntity::TET)
             {
                   for (int i=0;i<4;i++)
                   {
                           mVertex *V=(mVertex *)meshEntity->get(0,i);
                           fprintf(f,"%8d",(V->getId()+1));
                   }
                   fprintf(f,"\n");
             }
           }
         }
        }
       if (NBPRISM!=0)
       {
       fprintf(f,"penta6 \n%8d\n",NBPRISM);
       for(mMesh::iter it = theMesh->begin(3);it != theMesh->end(3); ++it)
         {
           AOMD::mEntity *meshEntity=*(it);
           mEntity::mType TYPELEM=meshEntity->getType ();
           AOMD::gEntity *geoEntity= meshEntity->getClassification();
           if (geoEntity==thegEntity)
           {
            if(TYPELEM == mEntity::PRISM)
             {
                   for (int i=0;i<6;i++)
                   {
                           mVertex *V=(mVertex *)meshEntity->get(0,i);
                           fprintf(f,"%8d",(V->getId()+1));
                   }
                   fprintf(f,"\n");
             }
           }
         }
        }
       if (NBHEX!=0)
       {
       fprintf(f,"hexa8 \n%8d\n",NBHEX);
       for(mMesh::iter it = theMesh->begin(3);it != theMesh->end(3); ++it)
         {
           AOMD::mEntity *meshEntity=*(it);
           mEntity::mType TYPELEM=meshEntity->getType ();
           AOMD::gEntity *geoEntity= meshEntity->getClassification();
           if (geoEntity==thegEntity)
           {
            if(TYPELEM == mEntity::HEX)
             {
                   for (int i=0;i<6;i++)
                   {
                           mVertex *V=(mVertex *)meshEntity->get(0,i);
                           fprintf(f,"%8d",(V->getId()+1));
                   }
                   fprintf(f,"\n");
             }
           }
         }
        }
     }
     if(f)fclose(f);
     */
}

mEntity *readGmshELM(istream &f, mMesh *theMesh, int recur, int &nb)
{
   int iNbNod, iTyp, iGrp, iElm, iNbSub, id;
   mVertex *nod[100];
   f >> iElm >> iTyp >> iGrp >> iNbSub >> iNbNod;
   for (int i = 0; i < iNbNod; i++)
   {
      f >> id;
      nod[i] = theMesh->getVertex(id);
      if (!nod[i])
      {
         char text[256];
         sprintf(text, "unknown verex id's %d, iNbNod=%d", id, i);
         throw new mException(__LINE__, __FILE__, text);
      }
   }

   nb++;

   mEntity *theEntity = nullptr;

   switch (iTyp)
   {
      case 4:
         theEntity = theMesh->createTetWithVertices(nod[0], nod[1], nod[2], nod[3], theMesh->getGEntity(iGrp, 3));
         break;
      case 5:
         theEntity = (mEntity *)theMesh->createHexWithVertices(nod[0], nod[1], nod[2], nod[3], nod[4], nod[5], nod[6], nod[7],
                                                               theMesh->getGEntity(iGrp, 3));
         break;
      case 6:
         theEntity = (mEntity *)theMesh->createPrismWithVertices(nod[0], nod[1], nod[2], nod[3], nod[4], nod[5],
                                                                 theMesh->getGEntity(iGrp, 3));
         break;
      case 2:
         theEntity = theMesh->createFaceWithVertices(nod[0], nod[1], nod[2], theMesh->getGEntity(iGrp, 2));
         break;
      case 3:
         theEntity = theMesh->createFaceWithVertices(nod[0], nod[1], nod[2], nod[3], theMesh->getGEntity(iGrp, 2));
         break;
      case 1:
         theEntity = theMesh->createEdge(nod[0], nod[1], theMesh->getGEntity(iGrp, 1));
         break;
      case 15:
         mVertex *v = nod[0];
         //	if(iGrp > 0)
         {
            v->classify(theMesh->getGEntity(iGrp, 0));
         }
         theEntity = v;
         break;
   }

      /// KH : read statically partitioned files


   {
      for (int i = 0; i < -iNbSub; i++)
      {
         mEntity *sub = readGmshELM(f, theMesh, recur + 1, nb);
         theEntity->add(sub);
         sub->setParent(theEntity);
      }
   }
   return theEntity;
}

// the same as above but for gmsh V2 msh file
mEntity *readGmshV2Elements(istream &f, mMesh *theMesh, int recur, int &nb)
{
   int iNbNod, iTyp, iGrp, iElm, iNbSub, id, iParDom;

   mVertex *nod[100];
   f >> iElm >> iTyp >> iParDom >> iGrp >> iNbSub;

   //   read nodes as long as there are nodes
   iNbNod = 0;
   string _line;
   getline(f, _line);
   vector<string> _words;
   _words.clear();
   char *pch;
   char _cline[_line.length() + 1];
   strcpy(_cline, _line.c_str());
   pch = strtok(_cline, " ");
   while (pch != nullptr)
   {
      _words.emplace_back(pch);
      pch = strtok(nullptr, " ");
   }
   iNbNod = _words.size();
   //   now we know how much nodes the element has

   for (int i = 0; i < iNbNod; i++)
   {
      id = atoi(_words[i].c_str());
      nod[i] = theMesh->getVertex(id);
      if (!nod[i])
      {
         char text[256];
         sprintf(text, "unknown vetrex id's %d, iNbNod=%d", id, i);
         throw new mException(__LINE__, __FILE__, text);
      }
   }

   nb++;

   mEntity *theEntity = nullptr;

   switch (iTyp)
   {
      case 4:
         theEntity = theMesh->createTetWithVertices(nod[0], nod[1], nod[2], nod[3], theMesh->getGEntity(iGrp, 3));
         break;
      case 5:
         theEntity = (mEntity *)theMesh->createHexWithVertices(nod[0], nod[1], nod[2], nod[3], nod[4], nod[5], nod[6], nod[7],
                                                               theMesh->getGEntity(iGrp, 3));
         break;
      case 6:
         theEntity = (mEntity *)theMesh->createPrismWithVertices(nod[0], nod[1], nod[2], nod[3], nod[4], nod[5],
                                                                 theMesh->getGEntity(iGrp, 3));
         break;
      case 2:
         theEntity = theMesh->createFaceWithVertices(nod[0], nod[1], nod[2], theMesh->getGEntity(iGrp, 2));
         break;
      case 3:
         theEntity = theMesh->createFaceWithVertices(nod[0], nod[1], nod[2], nod[3], theMesh->getGEntity(iGrp, 2));
         break;
      case 1:
         theEntity = theMesh->createEdge(nod[0], nod[1], theMesh->getGEntity(iGrp, 1));
         break;
      case 15:
         mVertex *v = nod[0];
         //	if(iGrp > 0)
         {
            v->classify(theMesh->getGEntity(iGrp, 0));
         }
         theEntity = v;
         break;
   }

      /// KH : read statically partitioned files


   {
      for (int i = 0; i < -iNbSub; i++)
      {
         mEntity *sub = readGmshELM(f, theMesh, recur + 1, nb);
         theEntity->add(sub);
         sub->setParent(theEntity);
      }
   }
   return theEntity;
}

/*
  A vertex can have several other vertices that are the same
  but in a different location. In AOMD, being the same means
  having the same iD. So, we can associate pairs of vertices
  that are supposed to be the same but with different
  locations, like mirror copies.

  associations like

  1 2
  3 4
  1 3
  2 4
  5 6

  are allowed.

  following algorithm will create two mMirrorVertex

  one containing vertex 1,2,3,4
  the other with vertices 5,6
*/

static void addAss(int source, std::set<int> &visited, std::set<int> &group, std::multimap<int, int> &vertexAssociations)
{
   if (group.find(source) == group.end())
   {
      group.insert(source);
      visited.insert(source);
      for (std::multimap<int, int>::iterator it = vertexAssociations.lower_bound(source);
           it != vertexAssociations.upper_bound(source); ++it)
      {
         addAss((*it).second, visited, group, vertexAssociations);
      }
   }
}

// nico
// static
void AOMD_Util::resolveAssociations(mMesh *theMesh, std::multimap<int, int> &vertexAssociations)
{
   std::set<int> visited;
   for (std::multimap<int, int>::iterator it = vertexAssociations.begin(); it != vertexAssociations.end(); ++it)
   {
      int source = (*it).first;
      std::set<int> group;
      if (visited.find(source) == visited.end())
      {
         addAss(source, visited, group, vertexAssociations);
      }
      if (group.size())
      {
         // nico	  mMirrorVertex *me = theMesh->createMirrorVertex();
         mMirrorVertex *me = new mMirrorVertex;
         for (std::set<int>::iterator it2 = group.begin(); it2 != group.end(); ++it2)
         {
            mVertex *v = theMesh->getVertex(*it2);
            if (!v)
            {
               char text[256];
               sprintf(text, "unknown verex id's %d in mirror vertex creation\n", *it2);
               throw new mException(__LINE__, __FILE__, text);
            }
            me->addCopy(v);
         }
      }
   }
}

void AOMD_Util::importDGFile(istream &f, mMesh *theMesh)
{
   char line[256];
   std::multimap<int, int> vertexAssociations;
   double gmshFormat;
   while (1)
   {
      do
      {
         f.getline(line, 256);
         if (f.eof()) break;
      } while (line[0] != '$');
      if (f.eof()) break;

      if (!strncmp(&line[1], "MeshFormat", 10))
      {
         int i1, i2;
         f >> gmshFormat >> i1 >> i2;
         //	cout << "gmsh mesh format " << gmshFormat << endl;
      }

      if ((!strncmp(&line[1], "NO", 2)) || (!strncmp(&line[1], "Nodes", 5)))
      {
         int NbNod;
         f >> NbNod;
         //	cout << "n nodes " <<  NbNod << endl;
         //	  printf("%d nodes\n",NbNod);
         for (int i = 0; i < NbNod; i++)
         {
            int iNod;
            double x, y, z;
            f >> iNod >> x >> y >> z;
            theMesh->createVertex(iNod, x, y, z, nullptr);
         }
      }
      if (!strncmp(&line[1], "SYM", 3))
      {
         int NbAss;
         f >> NbAss;
         for (int i = 0; i < NbAss; i++)
         {
            int v1, v2, nv;
            f >> nv >> v1;
            for (int j = 1; j < nv; j++)
            {
               f >> v2;
               std::pair<int const, int> pr(v1, v2);
               vertexAssociations.insert(pr);
            }
         }
         resolveAssociations(theMesh, vertexAssociations);
      }
      if (!strncmp(&line[1], "ELM", 3))
      {
         int NbElm;
         f >> NbElm;
         //	  printf("%d elms\n",NbElm);
         for (int i = 0; i < NbElm;)
         {
            //	      printf("%d | %d\n",i,NbElm);
            readGmshELM(f, theMesh, 0, i);
         }
      }
      if ((!strncmp(&line[1], "Elements", 8)) && (gmshFormat > 2))
      {
         int NbElm;
         f >> NbElm;
         //	cout << "n elem " <<  NbElm << endl;
         //	  printf("%d elms\n",NbElm);
         for (int i = 0; i < NbElm;)
         {
            //	      printf("%d | %d\n",i,NbElm);
            readGmshV2Elements(f, theMesh, 0, i);
         }
      }
      if (!strncmp(&line[1], "PRM", 3))
      {
         int NbParam, vert;
         double x, y, z;
         f >> NbParam;
         for (int i = 0; i < NbParam; i++)
         {
            f >> vert >> x >> y >> z;
            mVertex *vertex = theMesh->getVertex(vert);
            if (!vertex) throw;
            vertex->attachVector(getParametric(), Trellis_Util::mVector(x, y, z));
         }
      }
      if (!strncmp(&line[1], "FMO", 3))
      {
         int NbEmo, vert1, vert2, vert3, nbPt;
         double x, y, z;
         f >> NbEmo;
         // printf("There are total : %d  fmodes\n", NbEmo);
         for (int i = 0; i < NbEmo; i++)
         {
            f >> vert1 >> vert2 >> vert3 >> nbPt;
            mVertex *v_fmod_1 = theMesh->getVertex(vert1);
            mVertex *v_fmod_2 = theMesh->getVertex(vert2);
            mVertex *v_fmod_3 = theMesh->getVertex(vert3);
            if (!v_fmod_1 || !v_fmod_2 || !v_fmod_3)
            {
               char text[256];
               sprintf(text, "Unknown vertex %d %d %d <-> %p %p %p in face mode creation\n", vert1, vert2, vert3, v_fmod_1,
                       v_fmod_2, v_fmod_3);
               throw new mException(__LINE__, __FILE__, text);
            }

            mFace *face = theMesh->getTri(v_fmod_1, v_fmod_2, v_fmod_3);
            if (!face)
            {
               face = theMesh->createFaceWithVertices(vert1, vert2, vert3, theMesh->getGEntity(1000, 2));
            }
            // we move forward to pass qartic
            /* the face mode order is very important */
            if (nbPt > 1)
            {
               mAttachablePointVector *ap;
               ap = new mAttachablePointVector;
               double(*pnts)[3];
               pnts = new double[nbPt][3];
               int *index;
               index = new int[nbPt];
               for (int j = 0; j < 3; j++)
               {
                  mEntity *vt = face->get(0, j);
                  if (vt->equal(v_fmod_1))
                     index[0] = j;
                  else if (vt->equal(v_fmod_2))
                     index[1] = j;
                  else if (vt->equal(v_fmod_3))
                     index[2] = j;
               }
               {
                  for (int j = 0; j < 3; j++)
                  {
                     f >> pnts[index[j]][0] >> pnts[index[j]][1] >> pnts[index[j]][2];
                  }
               }
               {
                  for (int j = 0; j < 3; j++)
                  {
                     Trellis_Util::mPoint p = Trellis_Util::mPoint(pnts[j][0], pnts[j][1], pnts[j][2]);
                     (ap->v).push_back(p);
                  }
               }
               face->attachData(_ATTFMOD, ap);
               delete[] index;
               delete[] pnts;
            }
            else if (nbPt == 1)
            {
               mAttachablePointVector *ap;
               ap = new mAttachablePointVector;
               for (int j = 0; j < nbPt; j++)
               {
                  f >> x >> y >> z;
                  Trellis_Util::mPoint p = Trellis_Util::mPoint(x, y, z);
                  (ap->v).push_back(p);
               }
               face->attachData(_ATTFMOD, ap);
            }
         }
      }
      if (!strncmp(&line[1], "EMO", 3))
      {
         int NbEmo, vert1, vert2, nbPt;
         double x, y, z;
         f >> NbEmo;
         for (int i = 0; i < NbEmo; i++)
         {
            f >> vert1 >> vert2 >> nbPt;

            mVertex *v_emod_1 = theMesh->getVertex(vert1);
            mVertex *v_emod_2 = theMesh->getVertex(vert2);
            if (!v_emod_1 || !v_emod_2)
            {
               char text[256];
               sprintf(text, "Unknown vertex %d %d <-> %p %p in edge mode creation\n", vert1, vert2, v_emod_1, v_emod_2);
               throw new mException(__LINE__, __FILE__, text);
            }

            mEdge *e = theMesh->getEdge(v_emod_1, v_emod_2);
            if (!e)
            {
               e = theMesh->createEdge(vert1, vert2, theMesh->getGEntity(1000, 2));
            }
            // we move forward to pass qartic
            // assert (nbPt == 1);
            int sign = 0;
            if ((e->get(0, 0))->equal(v_emod_1)) sign = 1;
            mAttachablePointVector *ap;
            ap = new mAttachablePointVector;
            for (int j = 0; j < nbPt; j++)
            {
               f >> x >> y >> z;
               Trellis_Util::mPoint p = Trellis_Util::mPoint(x, y, z);
               if (sign == 1)
                  (ap->v).push_back(p);
               else if (sign == 0)
                  (ap->v).insert((ap->v).begin(), p);
            }
            printf("The mAOMD size of edge mode is: %lu", (ap->v).size());
            e->attachData(_ATTEMOD, ap);
         }
      }

      if (!strncmp(&line[1], "PRL", 3))
      {
         int NbConnx;
         f >> NbConnx;
         for (int i = 0; i < NbConnx; i++)
         {
            int id, to;
            f >> id >> to;
            //#ifdef PARALLEL
            theMesh->theOwnerManager->addNeighbor(id, to);
            //#endif
         }
      }
      if (!strncmp(&line[1], "FIN", 3)) break;
      do
      {
         f.getline(line, 256);
         if (f.eof()) break;
      } while (line[0] != '$');
   }
   classifyUnclassifiedVerices(theMesh);
}

// this function will not be valid any more since we unifies id and tag
// to be consistent with model
void AOMD_Util::setMappingFile(const char *fileName, mMesh *theMesh)
{
   FILE *f = fopen(fileName, "r");

   if (!f)
   {
      char text[256];
      sprintf(text, "Cannot open ids_2_tag Mapping file %s\n", fileName);
      throw mException(__LINE__, __FILE__, text);
   }
   while (!feof(f))
   {
      char typtag[12], typid[12], ctag[12], cid[12];
      int tag, id;
      fscanf(f, "%s %s %d %s %s %d", typtag, ctag, &tag, typid, cid, &id);
      if (!strcmp(typtag, "vertex"))
      {
         //	  printf("%d %d %d\n",id,0,tag);
         theMesh->getGEntity(tag, 0);
      }
      else if (!strcmp(typtag, "edge"))
      {
         // printf("%d %d %d\n",id,1,tag);
         theMesh->getGEntity(tag, 1);
      }
      else if (!strcmp(typtag, "face"))
      {
         // printf("%d %d %d\n",id,2,tag);
         theMesh->getGEntity(tag, 2);
      }
      else if (!strcmp(typtag, "region"))
      {
         // printf("%d %d %d\n",id,3,tag);
         theMesh->getGEntity(tag, 3);
      }
   }
   /*
     for (mMesh::giter it = theMesh->begin();
     it != theMesh->end(); ++ it)
     {
     printf("%d %d %d\n",(*it)->getId(),(*it)->getLevel(),(*it)->getTag());
     }
   */

   fclose(f);
}

void AOMD_Util::importVTUFile(const char *fName, mMesh *theMesh)
{
   cout << "importVTUFile : Opening " << fName << endl;
   std::ifstream ifs(fName);
   importVTUFile(ifs, theMesh);
   ifs.close();
}

void AOMD_Util::importVTUFile(ifstream &f, mMesh *theMesh)
{
   // Assumptions : only simplices
   // and only one type of element (linear)
   //

   const bool debug{false};
   const char *line_c;
   string line;
   int LocalistionNodes{0}, LocalistionConnec{0};
   unsigned long Position, PositionMaillage;

   unsigned int NbNodes, NbElements, DimX, DimT, DimXFictif;

   if (!f.is_open())
   {
      cout << "ERROR Unable to open vtu mesh file\n";
      exit(1);
   }

   // First, parse in order to get the position of the different parts of the file...
   // Indeed, there no specific order for the nodes and element parts with VTU format...
   // Inspired by ICITech code...
   while (!(f.eof()))
   {
      Position = f.tellg();
      getline(f, line);

      line_c = line.c_str();

      // Extract number of nodes and elements : NbNodes, NbElements
      if (line.find("NumberOfPoints") != string::npos)
      {
         // Extraction Nb noeuds et Nb elements
         long a, b;
         a = line.find("\"");
         b = line.find("\"", a + 1);
         sscanf(&(line_c[a + 1]), "%d", &NbNodes);  // ch_val.value_type

         a = line.find("\"", b + 1);
         sscanf(&(line_c[a + 1]), "%d", &NbElements);

         PositionMaillage = Position;
      }

      // Extraction nodes spatial dimensions : DimXFictif (always 3D I think)
      if (line.find("NumberOfComponents") != string::npos && line.find("Name") == string::npos)
      {
         long a, b;
         b = line.find("NumberOfComponents");
         a = line.find("\"", b + 1);
         sscanf(&(line_c[a + 1]), "%d", &DimXFictif);  // ch_val.value_type
      }

      // Extraction elements dimension : DimT
      // Then, deduce the min number of components for the nodes (2 in 2D, 3 in 3D) : DimX
      if (line.find("offsets") != string::npos)
      {
         Position = f.tellg();
         getline(f, line);
         line_c = line.c_str();

         sscanf(line_c, "%d", &DimT);
         DimX = DimT - 1;  // des simplex !!!
         // Gestion des fichiers de points (NbNodes=NbElements)
         if (DimT == 1) DimT = 4;
      }

      //<Greg
      // Extract the position of the section defining the nodes in the file
      if (line.find("<Points>") != string::npos)
      {
         getline(f, line);  // skip the next XML line
         LocalistionNodes = f.tellg();
      }

      // Extract the position of the section defining the connectivity in the file
      if (line.find("connectivity") != string::npos)
      {
         LocalistionConnec = f.tellg();
      }
      // Greg>
   }

   // Then, use the positions in order to fill the mesh
   // First, consider nodes
   f.clear();  // Whooooooo tricky ! If eof was raised, seekg can't go enywhere. We must clear f before or re-open it...
   f.seekg(LocalistionNodes);

   double xyzNode[3] = {0., 0., 0.};
   for (int iNode = 0; iNode < NbNodes; ++iNode)
   {
      for (int j = 0; j < DimX; ++j)
      {
         f >> xyzNode[j];
      }

      f.ignore('\n');  // Jump to next line

      if (debug) printf("Node %d : %f %f %f\n", iNode, xyzNode[0], xyzNode[1], xyzNode[2]);
      theMesh->createVertex(iNode, xyzNode[0], xyzNode[1], xyzNode[2], nullptr);
   }

   // Go to the connectivity...
   f.clear();  // se above
   f.seekg(LocalistionConnec);

   mFace *face;
   mTet *tet;
   unsigned int surfTag{1}, volTag{1};
   int connec[DimT];
   for (int iElt = 0; iElt < NbElements; ++iElt)
   {
      if (debug) cout << iElt << "connec ";

      for (int j = 0; j < DimT; j++)
      {
         f >> connec[j];
         if (debug) cout << connec[j] << " ";
      }

      if (debug) cout << endl;

      if (DimT == 3)
         face = theMesh->createFaceWithVertices(connec[0], connec[1], connec[2], theMesh->getGEntity(surfTag, 2));
      else
         tet = theMesh->createTetWithVertices(connec[0], connec[1], connec[2], connec[3], theMesh->getGEntity(volTag, 3));
   }

   f.close();
   classifyUnclassifiedVerices(theMesh);
}

void AOMD_Util::exportVTUFile(const char *fileName, const mMesh *theMesh, bool reduce_to_minimum)
{
   // We assume here that the nodes are numbered constiguously
   // Todo : add classified entities, add the physical tags as additional fields

   ofstream file(fileName);
   string Tab("  ");
   string Tab2 = Tab + Tab;
   string Tab3 = Tab + Tab + Tab;
   string Tab4 = Tab + Tab + Tab + Tab;
   string Tab5 = Tab + Tab + Tab + Tab + Tab;
   string Tab6 = Tab + Tab + Tab + Tab + Tab + Tab;
   string Tab7 = Tab + Tab + Tab + Tab + Tab + Tab + Tab;

   // Header
   file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
   file << Tab + "<UnstructuredGrid>" << endl;
   file << Tab2 + "<Piece NumberOfPoints=\"" << theMesh->size(0) << "\" NumberOfCells=\"" << theMesh->size(theMesh->getDim())
        << "\">" << endl;

   // Nodes
   file << Tab3 + "<Points>" << endl;
   file << Tab4 + "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

   // Check if node numbering begins at zero or one
   int nodesOffset{0};
   if (!theMesh->getVertex(0)) nodesOffset = 1;  //

   // can't use mesh iterators as they won't be ordered
   for (int i = 0 + nodesOffset; i < theMesh->size(0) + nodesOffset; ++i)
   {
      mVertex *v = theMesh->getVertex(i);
      file << Tab5 << v->point()(0) << " " << v->point()(1) << " " << v->point()(2) << endl;
   }

   file << Tab4 + "</DataArray>" << endl;
   file << Tab3 + "</Points>" << endl;

   // Elements
   // No problem with elements ordering here...
   file << Tab3 + "<Cells>" << endl;
   file << Tab4 + "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;

   unsigned int nbElts{0};
   if (!theMesh->size(3))
   {  // 2D ! --> Only linear triangles !!!
      nbElts = theMesh->size(2);
      for (mMesh::iter it = theMesh->begin(2); it != theMesh->end(2); ++it)
      {
         mEntity *e = *it;
         file << Tab5 << e->get(0, 0)->getId() - nodesOffset << " " << e->get(0, 1)->getId() - nodesOffset << " "
              << e->get(0, 2)->getId() - nodesOffset << endl;
      }
   }
   else
   {  // 3D ! --> Only linear tetrahedra !!!
      nbElts = theMesh->size(3);
      for (mMesh::iter it = theMesh->begin(3); it != theMesh->end(3); ++it)
      {
         mEntity *e = *it;
         file << Tab5 << e->get(0, 0)->getId() - nodesOffset << " " << e->get(0, 1)->getId() - nodesOffset << " "
              << e->get(0, 2)->getId() - nodesOffset << " " << e->get(0, 3)->getId() - nodesOffset << endl;
      }
   }

   file << Tab4 + "</DataArray>" << endl;
   file << Tab4 + "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;

   int eltOffset = 0;
   int cSize = theMesh->getDim() == 2 ? 3 : 4;
   file << Tab5;
   for (unsigned int i = 0; i < nbElts; ++i)
   {
      eltOffset += cSize;
      file << eltOffset << " ";
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
   file << Tab4 + "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;

   int eType = theMesh->getDim() == 2 ? 5 : 10;
   file << Tab5;
   for (unsigned int i = 0; i < nbElts; ++i)
   {
      file << eType << " ";
   }
   file << endl;
   file << Tab4 + "</DataArray>" << endl;
   file << Tab3 + "</Cells>" << endl;

   // Footer
   file << Tab2 + "</Piece>\n";
   file << Tab + "</UnstructuredGrid>\n";
   file << "</VTKFile>" << endl;
   file.close();
}

void AOMD_Util::exportICIFile(const char *fileName, const mMesh *theMesh, bool reduce_to_minimum)
{
   // We assume here that the nodes are numbered constiguously

   ofstream file(fileName);

   // Header
   int cSize = theMesh->getDim() == 2 ? 3 : 4;
   file << theMesh->size(0) << " " << theMesh->getDim() << " " << theMesh->size(theMesh->getDim()) << " " << cSize << endl;

   // Nodes

   // Check if node numbering begins at zero or one
   int nodesOffset{0};
   if (!theMesh->getVertex(0)) nodesOffset = 1;  //

   // can't use mesh iterators as they won't be ordered
   for (int i = 0 + nodesOffset; i < theMesh->size(0) + nodesOffset; ++i)
   {
      mVertex *v = theMesh->getVertex(i);
      file << v->point()(0) << " " << v->point()(1);
      if (theMesh->getDim() == 3) file << " " << v->point()(2);
      file << endl;
   }

   unsigned int nbElts{0};
   if (!theMesh->size(3))
   {  // 2D ! --> Only linear triangles !!!
      nbElts = theMesh->size(2);
      for (mMesh::iter it = theMesh->begin(2); it != theMesh->end(2); ++it)
      {
         mEntity *e = *it;
         file << e->get(0, 0)->getId() - nodesOffset << " " << e->get(0, 1)->getId() - nodesOffset << " "
              << e->get(0, 2)->getId() - nodesOffset << endl;
      }
   }
   else
   {  // 3D ! --> Only linear tetrahedra !!!
      nbElts = theMesh->size(3);
      for (mMesh::iter it = theMesh->begin(3); it != theMesh->end(3); ++it)
      {
         mEntity *e = *it;
         file << e->get(0, 0)->getId() - nodesOffset << " " << e->get(0, 1)->getId() - nodesOffset << " "
              << e->get(0, 2)->getId() - nodesOffset << " " << e->get(0, 3)->getId() - nodesOffset << endl;
      }
   }

   file.close();
}

}  // namespace AOMD
