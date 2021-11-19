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
#include "MeshModel.h"
#include "mEdge.h"
#include "mPoint.h"
#include "mMesh.h"
#include "mVertex.h"

#include "modeler.h"

using namespace AOMD;

void MeshModel_middlePoint (pMesh mesh, pEdge pE,     // an edge to be split
 double * pin)   // the boundary location of the vertex created in splitting
  
{

  mEdge *ed = (mEdge*) pE;
  Trellis_Util::mPoint p1 = ed->vertex(0)->point();
  Trellis_Util::mPoint p2 = ed->vertex(1)->point();
      
  pGEntity ge = pE->getClassification();   
  
//  AOMD::mMesh *mesh = ge->getMesh();  
  if (!mesh)throw;
  AOMD::mMesh *meshmodel = mesh->getMeshModel();

  Trellis_Util::mPoint pp ((p1+p2)*0.5);

  if (meshmodel)
    {

      mClassIterator it = meshmodel->begin (GEN_type(ge), GEN_tag(ge), GEN_type(ge));
      mClassIterator ite = meshmodel->end  (GEN_type(ge), GEN_tag(ge), GEN_type(ge));      
      double min_sum_d = 1.e20;
  /// model edge 
      if (GEN_type(ge) == 1)
	{	  
	  while (it != ite)
	    {
	      mEdge *e = (mEdge *) *it;
	      Trellis_Util::mPoint p1model = e->vertex(0)->point();
	      Trellis_Util::mPoint p2model = e->vertex(1)->point();
	      
	      Trellis_Util::mPoint pmid ((p1model+p2model)*0.5);

	      Trellis_Util::mVector v1 (p1, pmid);
	      Trellis_Util::mVector v2 (p2, pmid);

	      double sum = fabs(v1.normValue() - v2.normValue());

	      if (sum < min_sum_d)
		{
		  min_sum_d = sum;
		  pp = pmid;
		}
	      ++it;
	    }
	}
    }
  pin[0] = pp(0);
  pin[1] = pp(1);
  pin[2] = pp(2);
}
