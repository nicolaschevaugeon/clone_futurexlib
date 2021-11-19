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

/*<i> ****************************************************************************
 *
 *  AOMD/src/mDebugUtil.cc
 *  Created by Seegyoung Seol, on Mon Dec 08 2003, 10:26:37 EDT
 *
 *  File Content:debug utility function definitions
 *
 *************************************************************************** </i>*/

#include "mDebugUtil.h"
#include <cassert>
#include <iostream>
#include "AOMD_Internals.h"
#include "AOMD_cint.h"
#include "ParUtil.h"
#include "mAOMD.h"
#include "mEdge.h"
#include "mEntity.h"
#include "mException.h"
#include "mFace.h"
#include "mHex.h"
#include "mMesh.h"
#include "mPoint.h"
#include "mPrism.h"
#include "mTet.h"
#include "mVertex.h"



#ifndef SIM
#include "modeler.h"
#else
#include "SimModel.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
using std::cout;
using std::endl;
using std::list;
using std::pair;
using std::set;
using std::sort;
using std::string;
using std::vector;

namespace AOMD
{
#ifndef SIM
// ***********************************************************
bool verify_mesh(mMesh* mesh)
// ***********************************************************
{
   pmEntity* pe;
   mEntity* vertex;
   mEntity* edge;
   mEntity* face;
   mEntity* region;
   int gLevel;
   int manifold = 1;
   int numAdjE;
   int numAdjF;
   int numAdjR;
   int mypid = P_pid();
   int meshDim = M_globalMaxDim(mesh);
   int ret = 1; /* true */



   // check vertices
   for (mMesh::iterall it = mesh->beginall(0); it != mesh->endall(0); ++it)
   {
      vertex = *it;
      numAdjE = V_numEdges(vertex);
      if (meshDim == 3 && numAdjE < 3)
      {
         cout << "(" << mypid << ") AOMD WARNING: " << numAdjE << " adj. edges with " << vertex->getUid() << endl;
         vertex->print();
         if (numAdjE == 0) manifold = 0;
      }
      if (meshDim == 2 && numAdjE < 2)
      {
         cout << "(" << mypid << ") AOMD WARNING: no adj. edges with " << vertex->getUid() << endl;
         vertex->print();
         manifold = 0;
      }

      pe = vertex->getPClassification();

      assert(!pe);

   }

   // check edges
   for (mMesh::iterall it = mesh->beginall(1); it != mesh->endall(1); ++it)
   {
      edge = *it;
      numAdjF = E_numFaces(edge);

      if (meshDim == 3)
      {
         switch (numAdjF)
         {
            case 0:
               cout << "(" << mypid << ") AOMD WARNING: " << E_numFaces(edge) << " adj. faces with " << edge->getUid() << endl;
               edge->print();
               manifold = 0;
               break;
            case 1:
               if (!EN_duplicate(edge))
               {
                  cout << "(" << mypid << ") AOMD WARNING: " << E_numFaces(edge) << " adj. faces with " << edge->getUid() << endl;
                  edge->print();
                  manifold = 0;
               }
               break;
            default:
               break;
         }
         gLevel = GEN_type(edge->getClassification());
         if (gLevel == 1 || gLevel == 2)
         {
            int bf_counter = 0;  // number of boundary faces adjacent
            for (int i = 0; i < numAdjF; ++i)
            {
               face = edge->get(2, i);
               if (GEN_type(face->getClassification()) == 2) bf_counter++;
            }
            if (bf_counter != 2 && !EN_duplicate(edge))
            {
               edge->print();
               cout << "(" << mypid << ") AOMD WARNING: non-manifold model? # adj. bounding faces=" << bf_counter << endl;
            }
            if (!EN_duplicate(edge))
               assert(bf_counter == 2);
            else
            {
               //	  assert(bf_counter==1);
            }
         }
      }
      else if (meshDim == 2)  // edge in 2D
      {
         if (numAdjF == 0)
         {
            cout << "(" << mypid << ") AOMD WARNING: no adj. faces with " << edge->getUid() << endl;
            edge->print();
         }

         gLevel = GEN_type(edge->getClassification());
         switch (gLevel)
         {
            case 1:  // model edge
               if (numAdjF != 1)
               {
                  cout << "(" << mypid << ") AOMD WARNING: " << numAdjF << " faces with ";
                  edge->print();
                  manifold = 0;
               }
               break;
            case 2:
               if (edge->getPClassification())  // CB edge
                  assert(numAdjF == 1);
               else  // non-CB edge
                  assert(numAdjF == 2);
               break;
            default:
               cout << "(" << mypid << ") AOMD ERROR: mesh edge is wrong classified\n";
               cout << "Mesh Dim = " << meshDim << endl;
               edge->print();
               ret = 0;
         }  // switch
      }     // else

      pe = edge->getPClassification();

      assert(!pe);

   }  // edges

   // check faces
   for (mMesh::iterall it = mesh->beginall(2); it != mesh->endall(2); ++it)
   {
      face = *it;

      numAdjR = F_numRegions(face);

      if (meshDim == 3)
      {
         gLevel = GEN_type(face->getClassification());
         switch (gLevel)
         {
            case 2:  // face on model face
               if (numAdjR != 1)
               {
                  cout << "(" << mypid << ") # adj. regions with " << face->getUid() << ": " << numAdjR << endl;
                  manifold = 0;
               }
               break;
            case 3:                     // interior face
               if (EN_duplicate(face))  // CB face
               {
                  if (numAdjR != 1) face->print();
                  assert(numAdjR == 1);
               }
               else  // non-CB face
               {
                  if (F_numRegions(face) != 2) face->print();
                  assert(F_numRegions(face) == 2);
               }
               break;
            default:
               cout << "(" << mypid << ") AOMD ERROR: mesh face is wrong classified. ";
               face->print();
               ret = 0;
         }
      }  // face in 3D
      else if (meshDim == 2)
      {
         if (numAdjR != 0)
         {
            cout << "(" << mypid << ") AOMD ERROR: adj. regions with " << face->getUid() << " in 2D mesh\n";
            ret = 0;
         }
         assert(F_numRegions(face) == 0);
         assert(GEN_type(face->getClassification()) == 2);
      }  // face in 2D

      pe = face->getPClassification();

      assert(!pe);

   }

   // check regions
   for (mMesh::iterall it = mesh->beginall(3); it != mesh->endall(3); ++it)
   {
      region = *it;
      assert(!(region->getPClassification()));
      assert(GEN_type(region->getClassification()) == 3);
#ifdef DEBUG
      if (R_Volume(region) < 0.0)
      {
         cout << "AOMD FATAL: Negative volume region " << region->getUid() << " found!\n   => ";
         region->print();
      }
#endif
      assert(R_Volume(region) > 0.0);
   }


   if (P_getMinInt(manifold) == 0) ParUtil::Instance()->Msg(ParUtil::INFO, "\n* Is Your Model Non-Manifold?\n");

   int globalRetVal;
   globalRetVal = P_getMinInt(ret);
   if (globalRetVal == 1)
      return true;
   else
      return false;
}  // verify_mesh

int mMesh::numUniqueEntities(int type)
{
   if (ParUtil::Instance()->size() == 1) return size(type);
   int mypid = ParUtil::Instance()->rank();
   int numEnt = 0;
   for (iterall it = beginall(type); it != endall(type); ++it)
      if ((*it)->getOwner() == mypid) numEnt++;
   return numEnt;
}

void mMesh::printAll() const
{
   for (int i = 0; i <= getDim(); ++i) print(i);
}

void mMesh::print(int dim) const
{


   for (iterall it = beginall(dim); it != endall(dim); ++it) (*it)->print();
}

void mMesh::printNumEntities()
{
   for (int i = 0; i <= 3; ++i)
   {
      vector<int> output;
      M_numEntitiesOwned(this, i, output);
      if (P_pid() == 0)
      {
         int sum = 0;
         for (int j = 0; j < P_size(); ++j) sum += output[j];
         cout << "# dim[" << i << "] = " << sum << " (";
         for (int j = 0; j < P_size(); ++j)
         {
            cout << output[j];
            if (j != P_size() - 1) cout << ",";
         }
         cout << ")" << endl;
      }
   }
}

void mMesh::printNumEntities_NonCfm()
{

   cout << "\n***  M_PrintNumEntities  ***\n";
   cout << "   numVt = " << numUniqueEntities(0) << endl;
   cout << "   numEg = " << numUniqueEntities(1) << endl;
   cout << "   numFc = " << numUniqueEntities(2) << endl;
   cout << "   numRg = " << numUniqueEntities(3) << endl;
   cout << "*****************************\n";

}

#endif /* SIM */

void mVertex::print() const
{

   cout << "en_info: " << getUid() << "[" << p(0) << "," << p(1) << "," << p(2) << "], "
        << "(iD " << iD << ", gLevel " << GEN_type(theClassification) << ", gTag " << GEN_tag(theClassification) << ", gId "
        << GEN_id(theClassification) << ")\n";

}

void mEdge::print() const
{

   cout << "en_info: " << getUid() << "(iD " << iD << ", gLevel " << GEN_type(theClassification) << ", gTag "
        << GEN_tag(theClassification) << ")\n";

}

void mFace::print() const
{

   cout << "en_info: " << getUid() << "(iD " << iD << ", gLevel " << GEN_type(theClassification) << ", gTag "
        << GEN_tag(theClassification) << ")\n";

}

void mTet::print() const
{

   cout << "en_info: " << getUid() << "(iD " << iD << ", gLevel " << GEN_type(theClassification) << ", gTag "
        << GEN_tag(theClassification) << ")\n";

}

void mHex::print() const
{

   cout << "en_info: " << getUid() << "(iD " << iD << ", gLevel " << GEN_type(theClassification) << ", gTag "
        << GEN_tag(theClassification) << ")\n";

}

void mPrism::print() const
{

   cout << "en_info: " << getUid() << "(iD " << iD << ", gLevel " << GEN_type(theClassification) << ", gTag "
        << GEN_tag(theClassification) << ")\n";

}
}  // namespace AOMD
