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

#ifndef _MESH_MODEL_
#define _MESH_MODEL_

#ifdef SIM
#include "MeshSim.h" 
#include "MeshSimInternal.h"
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void MeshModel_middlePoint
  ( pMesh,
  pEdge pE,       // an edge to be split
    double * pin );   // the boundary location of the vertex created in splitting
#ifdef __cplusplus
}
#endif
#endif
