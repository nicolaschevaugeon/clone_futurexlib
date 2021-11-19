/********************************************************************************

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

*********************************************************************************/
#include <cassert>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "AOMD_Internals.h"
#include "ParUtil.h"
#include "mAOMD.h"
#include "mEdge.h"
#include "mException.h"
#include "mFace.h"
#include "mMesh.h"
#include "mTet.h"
#include "mVertex.h"

#include "AOMD.h"
#include "SGModel.h"
#include "modeler.h"
using std::cout;
using std::endl;
using std::istream;
using std::ofstream;
using std::ostream;

namespace AOMD
{
void classifyUnclassifiedVerices(mMesh *m);

// **********************************************************
void AOMD_Util::importSmsFile(const char *fName, mMesh *theMesh)
// **********************************************************
{
#ifndef SIM
   FILE *in = fopen(fName, "r");
   if (!in) throw new mException(__LINE__, __FILE__, "unable to open file in loadSmsFile");

   char line[256];
   int i;

   int NbRegions, NbFaces, NbEdges, NbVertices, NbPoints, GEntityType, GEntityId, EntityNbConnections, Dummy, Edge1, Edge2, Edge3,
       Edge4, Face[6], nbPts;
   int VertexId1, VertexId2, NbEdgesOnFace, NbFacesOnRegion;
   double x, y, z;

   fscanf(in, "%s %d", line, &Dummy);
   typedef GEntity *(*entityBy_FP)(SGModel *, int, int);
   entityBy_FP fp = GM_entityByTag;

   if (Dummy == 2)
   {
      fp = GM_entityByID;
      theMesh->setGEntity_FP(GM_entityByTag);
   }

   fscanf(in, "%d %d %d %d %d", &NbRegions, &NbFaces, &NbEdges, &NbVertices, &NbPoints);
   //    printf("reading %d vertices ...\n",NbVertices);
   double u, v;
   int patch;

   for (i = 0; i < NbVertices; i++)
   {
      fscanf(in, "%d", &GEntityId);
      if (GEntityId)
      {
         fscanf(in, "%d %d %lf %lf %lf", &GEntityType, &EntityNbConnections, &x, &y, &z);
         mVertex *vv = theMesh->createVertex(i + 1, x, y, z, theMesh->getGEntity(GEntityId - 1, GEntityType, fp));

         switch (GEntityType)
         {
            case 0:
               break;
            case 1:
               fscanf(in, "%le", &u);
               vv->attachVector(getParametric(), Trellis_Util::mVector(u, 0, 0));
               break;
            case 2:
               fscanf(in, "%le %le %d", &u, &v, &patch);
               vv->attachVector(getParametric(), Trellis_Util::mVector(u, v, patch));
               break;
            case 3:
               break;
            default:
               if (!in) throw new mException(__LINE__, __FILE__, "unknown tag in sms loadSmsFile");
         }
      }
   }

   int *edges = new int[2 * NbEdges];

   //    printf("reading %d edges ...\n",NbEdges);
   for (i = 0; i < NbEdges; i++)
   {
      fscanf(in, "%d", &GEntityId);
      if (GEntityId)
      {
         fscanf(in, "%d %d %d %d %d", &GEntityType, &VertexId1, &VertexId2, &EntityNbConnections, &nbPts);
         edges[2 * i] = VertexId1;
         edges[2 * i + 1] = VertexId2;

         if (GEntityType < 2) theMesh->createEdge(VertexId1, VertexId2, theMesh->getGEntity(GEntityId - 1, GEntityType, fp));

         for (int j = 0; j < nbPts; j++)
         {
            switch (GEntityType)
            {
               case 0:
                  break;
               case 1:
                  fscanf(in, "%le", &u);
                  break;
               case 2:
                  fscanf(in, "%le %le %d", &u, &v, &patch);
                  break;
               case 3:
                  break;
               default:
                  delete[] edges;
                  throw new mException(__LINE__, __FILE__, "unknown tag in sms loadSmsFile");
            }
         }
      }
   }

   int *faces = new int[4 * NbFaces];

   //    printf("reading %d faces ...\n",NbFaces);
   for (i = 0; i < NbFaces; i++)
   {
      fscanf(in, "%d", &GEntityId);
      if (GEntityId)
      {
         fscanf(in, "%d %d", &GEntityType, &NbEdgesOnFace);

         if (NbEdgesOnFace == 3)
         {
            fscanf(in, "%d %d %d %d", &Edge1, &Edge2, &Edge3, &nbPts);
            int iVertex1 = (Edge1 > 0) ? edges[2 * (abs(Edge1) - 1)] : edges[2 * (abs(Edge1) - 1) + 1];
            int iVertex2 = (Edge2 > 0) ? edges[2 * (abs(Edge2) - 1)] : edges[2 * (abs(Edge2) - 1) + 1];
            int iVertex3 = (Edge3 > 0) ? edges[2 * (abs(Edge3) - 1)] : edges[2 * (abs(Edge3) - 1) + 1];
            faces[4 * i] = iVertex1;
            faces[4 * i + 1] = iVertex2;
            faces[4 * i + 2] = iVertex3;

            if (GEntityType < 3)
               theMesh->createFaceWithVertices(theMesh->getVertex(iVertex1), theMesh->getVertex(iVertex2),
                                               theMesh->getVertex(iVertex3), theMesh->getGEntity(GEntityId - 1, GEntityType, fp));
         }
         else if (NbEdgesOnFace == 4)
         {
            fscanf(in, "%d %d %d %d %d", &Edge1, &Edge2, &Edge3, &Edge4, &nbPts);
            int iVertex1 = (Edge1 > 0) ? edges[2 * (abs(Edge1) - 1)] : edges[2 * (abs(Edge1) - 1) + 1];
            int iVertex2 = (Edge2 > 0) ? edges[2 * (abs(Edge2) - 1)] : edges[2 * (abs(Edge2) - 1) + 1];
            int iVertex3 = (Edge3 > 0) ? edges[2 * (abs(Edge3) - 1)] : edges[2 * (abs(Edge3) - 1) + 1];
            int iVertex4 = (Edge4 > 0) ? edges[2 * (abs(Edge4) - 1)] : edges[2 * (abs(Edge4) - 1) + 1];
            faces[4 * i] = iVertex1;
            faces[4 * i + 1] = iVertex2;
            faces[4 * i + 2] = iVertex3;
            faces[4 * i + 3] = iVertex4;

            if (GEntityType < 3)
               theMesh->createFaceWithVertices(theMesh->getVertex(iVertex1), theMesh->getVertex(iVertex2),
                                               theMesh->getVertex(iVertex3), theMesh->getVertex(iVertex4),
                                               theMesh->getGEntity(GEntityId - 1, GEntityType, fp));
         }
         else
         {
            delete[] faces;
            delete[] edges;
            throw new mException(__LINE__, __FILE__, "wrong number of edges on a face in loadSmsFile");
         }
         for (int j = 0; j < nbPts; j++)
         {
            switch (GEntityType)
            {
               case 0:
                  break;
               case 1:
                  fscanf(in, "%le", &u);
                  break;
               case 2:
                  fscanf(in, "%le %le %d", &u, &v, &patch);
                  break;
               case 3:
                  break;
               default:
                  delete[] faces;
                  delete[] edges;
                  throw new mException(__LINE__, __FILE__, "unknown tag in sms loadSmsFile");
            }
         }
      }
   }

   delete[] edges;

   //    printf("reading %d regions ...\n",NbRegions);
   for (i = 0; i < NbRegions; i++)
   {
      fscanf(in, "%d", &GEntityId);
      if (GEntityId)
      {
         fscanf(in, "%d", &NbFacesOnRegion);
         if (NbFacesOnRegion == 4)
         {
            fscanf(in, "%d %d %d %d %d", &Face[0], &Face[1], &Face[2], &Face[3], &Dummy);

            int iVertex1 = faces[4 * (abs(Face[0]) - 1)];
            int iVertex2 = faces[4 * (abs(Face[0]) - 1) + 1];
            int iVertex3 = faces[4 * (abs(Face[0]) - 1) + 2];
            int iVertex4;

            int iVertex = faces[4 * (abs(Face[1]) - 1)];
            if (iVertex != iVertex1 && iVertex != iVertex2 && iVertex != iVertex3)
            {
               iVertex4 = iVertex;
            }

            iVertex = faces[4 * (abs(Face[1]) - 1) + 1];
            if (iVertex != iVertex1 && iVertex != iVertex2 && iVertex != iVertex3)
            {
               iVertex4 = iVertex;
            }

            iVertex = faces[4 * (abs(Face[1]) - 1) + 2];
            if (iVertex != iVertex1 && iVertex != iVertex2 && iVertex != iVertex3)
            {
               iVertex4 = iVertex;
            }

            theMesh->createTetWithVertices(theMesh->getVertex(iVertex1), theMesh->getVertex(iVertex2),
                                           theMesh->getVertex(iVertex3), theMesh->getVertex(iVertex4),
                                           theMesh->getGEntity(GEntityId - 1, 3, fp));
         }
         if (NbFacesOnRegion == 6)
         {
            fscanf(in, "%d %d %d %d %d %d %d", &Face[0], &Face[1], &Face[2], &Face[3], &Face[4], &Face[5], &Dummy);

            int iVertex1 = faces[4 * (abs(Face[0]) - 1)];
            int iVertex2 = faces[4 * (abs(Face[0]) - 1) + 1];
            int iVertex3 = faces[4 * (abs(Face[0]) - 1) + 2];
            int iVertex4 = faces[4 * (abs(Face[0]) - 1) + 3];
            ;
            int iVertex5 = faces[4 * (abs(Face[5]) - 1)];
            int iVertex6 = faces[4 * (abs(Face[5]) - 1) + 1];
            int iVertex7 = faces[4 * (abs(Face[5]) - 1) + 2];
            int iVertex8 = faces[4 * (abs(Face[5]) - 1) + 3];

            theMesh->createHexWithVertices(
                theMesh->getVertex(iVertex1), theMesh->getVertex(iVertex2), theMesh->getVertex(iVertex3),
                theMesh->getVertex(iVertex4), theMesh->getVertex(iVertex5), theMesh->getVertex(iVertex6),
                theMesh->getVertex(iVertex7), theMesh->getVertex(iVertex8), theMesh->getGEntity(GEntityId - 1, 3, fp));
         }
      }
   }
   delete[] faces;
   fclose(in);
   classifyUnclassifiedVerices(theMesh);
#endif
}

// **********************************************************
void AOMD_Util::exportSmsFile(const char *fName, const mMesh *theMesh)
// **********************************************************
{
   // ---  this need to be modified later to be compatible with DMUM
   // --- E.S. Seol
   double t1 = ParUtil::Instance()->wTime();
   ofstream ofs(fName);
   std::map<int, int> numberid;
   double loc[3];
   pVertex vertex;
   pRegion region;
   pEdge edge;
   pFace face;

   int gentitytype, gentid, i, j, numEdges, numFaces, ID;
   pGEntity classifyon;

   int NbRegions, NbFaces, NbEdges, NbVertices, NbPoints;
   NbRegions = theMesh->size(3);
   NbFaces = theMesh->size(2);
   NbEdges = theMesh->size(1);
   NbVertices = theMesh->size(0);
   NbPoints = NbVertices;
   double low, high;
   /*
       cout<<"NbRegions="<<NbRegions<<endl;
       cout<<"NbFaces="<<NbFaces<<endl;
       cout<<"NbEdges="<<NbEdges<<endl;
       cout<<"NbVertices="<<NbVertices<<endl;
   */

   if (!strcmp((char *)(((SGModel *)theMesh->getSolidModel())->modeler()), "acis"))
   {
      pGModel model = theMesh->getSolidModel();

      //     if (ParUtil::Instance()->rank()==0)
      cout << "export mesh into " << fName << " version 4";
      ofs << "sms 4\n";  // old format
      // write out # model entities
      ofs << GM_numVertices(model) << " " << GM_numEdges(model) << " " << GM_numFaces(model) << " " << GM_numRegions(model)
          << endl;

      double xyz[3], xyz2[3];
      pGVertex gvertex;
      pGEdge gedge;
      pGFace gface;
      pGRegion gregion;
      pPList gevertices, gefaces;
      pPList gfedges, gfvertices;
      pPList grfaces;
      void *iter;
      // write out gvertices
      GVIter gvit = GM_vertexIter(model);
      while ((gvertex = GVIter_next(gvit)))
      {
         GV_point(gvertex, xyz);
         ofs << GEN_type((pGEntity)gvertex) << " " << GEN_tag((pGEntity)gvertex) << " " << xyz[0] << " " << xyz[1] << " "
             << xyz[2] << "\n";
      }

      GVIter_delete(gvit);

      // white out gedges
      GEIter geit = GM_edgeIter(model);
      while ((gedge = GEIter_next(geit)))
      {
         ofs << GEN_type((pGEntity)gedge) << " " << GEN_tag((pGEntity)gedge) << " ";

         // write evertices
         gevertices = GE_vertices(gedge);
         ofs << PList_size(gevertices) << " ";

         iter = nullptr;
         while ((gvertex = (pGVertex)PList_next(gevertices, &iter))) ofs << GEN_tag((pGEntity)gvertex) << " ";
         PList_delete(gevertices);

         gefaces = GE_faces(gedge);
         ofs << PList_size(gefaces) << " ";

         iter = nullptr;
         while ((gface = (pGFace)PList_next(gefaces, &iter)))
         {
            gfvertices = GF_vertices(gface);
            ofs << PList_size(gfvertices) << " ";

            void *iter2 = nullptr;
            while ((gvertex = (pGVertex)PList_next(gfvertices, &iter2))) ofs << GEN_tag((pGEntity)gvertex) << " ";
            PList_delete(gfvertices);
         }
         ofs << "\n";
         PList_delete(gefaces);
      }
      GEIter_delete(geit);

      // white out gfaces
      GFIter gfit = GM_faceIter(model);
      while ((gface = GFIter_next(gfit)))
      {
         ofs << GEN_type((pGEntity)gface) << " " << GEN_tag((pGEntity)gface) << " ";

         gfedges = GF_edges(gface);
         ofs << PList_size(gfedges) << " ";

         iter = nullptr;
         while ((gedge = (pGEdge)PList_next(gfedges, &iter))) ofs << GEN_tag((pGEntity)gedge) << " ";
         PList_delete(gfedges);

         gfvertices = GF_vertices(gface);
         ofs << PList_size(gfvertices) << " ";
         iter = nullptr;
         while ((gvertex = (pGVertex)PList_next(gfvertices, &iter))) ofs << GEN_tag((pGEntity)gvertex) << " ";
         GF_parRange(gface, 0, &low, &high);
         ofs << low << " " << high << " ";
         GF_parRange(gface, 1, &low, &high);
         ofs << low << " " << high << " ";

         ofs << "\n";
         PList_delete(gfvertices);
      }  // while
      GFIter_delete(gfit);

      // white out gregions
      GRIter grit = GM_regionIter(model);
      while ((gregion = GRIter_next(grit)))
      {
         ofs << GEN_type((pGEntity)gregion) << " " << GEN_tag((pGEntity)gregion) << " ";
         grfaces = GR_faces(gregion);
         ofs << PList_size(grfaces) << " ";

         iter = nullptr;
         while ((gface = (pGFace)PList_next(grfaces, &iter))) ofs << GEN_tag((pGEntity)gface) << " ";
         ofs << "\n";
         PList_delete(grfaces);
      }  // while
      GRIter_delete(grit);
   }  // if acis model
   else
   {
      //    if (ParUtil::Instance()->rank()==0)
      cout << "export mesh into " << fName << " version 1";
   }

   ofs << "sms 2\n";  // new format
   Trellis_Util::mVector vec;

   ofs << NbRegions << " " << NbFaces << " " << NbEdges << " " << NbVertices << " " << NbPoints << endl;

   if (NbVertices > 0)
   {
      i = 1;

      for (mMesh::iterall it = theMesh->beginall(0); it != theMesh->endall(0); ++it)
      {
         vertex = *it;
         numberid.insert(std::make_pair(EN_id(vertex), i++));
         V_coord(vertex, loc);
         classifyon = V_whatIn(vertex);
         gentitytype = GEN_type(classifyon);
         gentid = GEN_id(classifyon);

         ofs << gentid + 1 << " " << gentitytype << " " << V_numEdges(vertex) << endl
             << loc[0] << " " << loc[1] << " " << loc[2] << " ";
         if (vertex->getData(AOMD_Util::Instance()->getParametric()))
         {
            vec = vertex->getAttachedVector(AOMD_Util::Instance()->getParametric());
            switch (gentitytype)
            {
               case 1:
                  ofs << vec(0);
                  break;
               case 2:
                  ofs << vec(0) << " " << vec(1) << " " << vec(2);
                  break;
               default:
                  break;
            }
         }
         else
         {
            switch (gentitytype)
            {
               case 1:
                  ofs << " 0";
                  break;
               case 2:
                  ofs << " 0"
                      << " 0"
                      << " 0";
                  break;
               default:
                  break;
            }
         }
         ofs << endl;
      }  // for
   }     // if NbVertices

   if (NbEdges > 0)
   {
      i = 1;
      for (mMesh::iterall it = theMesh->beginall(1); it != theMesh->endall(1); ++it)
      {
         edge = *it;
         classifyon = E_whatIn(edge);
         gentitytype = GEN_type(classifyon);
         gentid = GEN_id(classifyon);
         ofs << gentid + 1 << " " << gentitytype << " ";
         EN_setID(edge, i++);
         for (j = 0; j < 2; j++)
         {
            vertex = E_vertex(edge, j);
            ID = EN_id(vertex);
            std::map<int, int>::const_iterator p = numberid.find(ID);
            ofs << (*p).second << " ";
         }
         ofs << E_numFaces(edge) << " "
             << "0";
         ofs << endl;
      }  // while
   }     // if NbEdges

   if (NbFaces > 0)
   {
      i = 1;
      for (mMesh::iterall it = theMesh->beginall(2); it != theMesh->endall(2); ++it)
      {
         face = *it;
         classifyon = F_whatIn(face);
         gentitytype = GEN_type(classifyon);
         gentid = GEN_id(classifyon);
         ofs << gentid + 1 << " " << gentitytype << " ";
         EN_setID(face, i++);
         numEdges = F_numEdges(face);
         ofs << numEdges << " ";
         for (j = 0; j < numEdges; j++)
         {
            edge = F_edge(face, j);
            ID = EN_id(edge);
            if (!F_edgeDir(face, j)) ID = -ID;
            ofs << ID << " ";
         }
         ofs << "0";  // NbPoints
         ofs << endl;
      }  // while
   }     // if NbFaces

   if (NbRegions > 0)
   {
      for (mMesh::iterall it = theMesh->beginall(3); it != theMesh->endall(3); ++it)
      {
         region = *it;
         classifyon = (pGEntity)R_whatIn(region);
         gentid = GEN_id(classifyon);
         numFaces = R_numFaces(region);
         ofs << gentid + 1 << " " << numFaces << " ";
         for (j = 0; j < numFaces; j++)
         {
            face = R_face(region, j);
            ID = EN_id(face);
            if (!R_faceDir(region, j)) ID = -ID;
            ofs << ID << " ";
         }
         ofs << "0";
         ofs << endl;
      }  // while
   }     // if

   ofs.close();
   double t2 = ParUtil::Instance()->wTime();
   //  if (ParUtil::Instance()->rank()==0)
   cout << " (" << t2 - t1 << " sec)\n";
}

}  // namespace AOMD
