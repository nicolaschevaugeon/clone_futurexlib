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
#include "AOMD.h"
#include "AOMDInternals.h"
#include "AOMD_cint.h"

using namespace AOMD;

pEdge E_exist(pVertex v1, pVertex v2)
{
  /*
  pVertex   v1    ; vertex  1 (I)
  pVertex   v2    ; vertex  2 (I)
  */

  pEdge edge ;
  pVertex evert[2] ;
  int nedge=0,i ;

  if (v1 == v2) return((pEdge)nullptr);

  nedge = V_numEdges(v1) ;

  if(nedge == 0) return((pEdge)nullptr);
 
  for(i=0;i<nedge;i++)
  {
    edge = V_edge(v1, i);
    /* check if both vertices of this edge match the given vertices */
    evert[0] = E_vertex(edge, 0);
    if(v2 == evert[0]) return(edge) ;
    evert[1] = E_vertex(edge, 1);
    if(v2 == evert[1]) return(edge) ;
  } 

  /* if here, then no valid connection found */
  return((pEdge)nullptr);
}

pFace F_exist(int type, pEntity e1, pEntity e2, pEntity e3, pEntity e4)
{
 pEntity ents[4] ;
 pFace face ;
 int count, i;
 pPList entities, faces;
 void *tmp = nullptr;

 ents[0] = e1;
 ents[1] = e2;
 ents[2] = e3;

 if (e4) {
   ents[3] = e4;
   count = 4;
 } else
   count = 3;

 /* if ents are vertices get the vertices of the face
 */

 int Tedge = 1;
 int Tvertex =0;

 if(type>Tedge)
   return((pFace)nullptr);

 if( type == Tvertex )
 {
   faces = V_faces((pVertex) ents[0]);
   /* For each face connected to ents[0] see if the other vertices are 
      connected to the face
   */
   
   for (;(face = (pFace)PList_next(faces, &tmp));) 
   {
     entities = F_vertices(face, 1);
     for( i = 1; i < count; i++)
       if (!PList_inList(entities, ents[i]))
	 break;
     PList_delete(entities);
     if (i == count){
       PList_delete(faces);/* Face is connected to all entities */
       return face;
     }
   }
 }
 else if (type == Tedge)
 {
   faces = E_faces((pEdge) ents[0]);
   /* For each face connected to ents[0] see if the other edges are 
      connected to the face
   */
   
   for (;(face = (pFace)PList_next(faces, &tmp));) {
     entities = F_edges(face, 1, (pVertex) nullptr);
     for( i = 1; i < count; i++)
       if (!PList_inList(entities, ents[i]))
	 break;
     PList_delete(entities);
     if (i == count) {
       PList_delete(faces);/* Face is connected to all entities */
       return face;
     }
   }
 }
 PList_delete(faces);

 return (pFace)nullptr ;
}

pRegion R_exist(int type, pEntity e1, pEntity e2, pEntity e3, 
                         pEntity e4, pEntity e5, pEntity e6)
{
 pEntity ents[6] ;
 pRegion region;
 
 int i, count;

 pPList entities, regions;
 void *tmp = nullptr;

 ents[0] = e1;
 ents[1] = e2;
 ents[2] = e3;
 ents[3] = e4;
 count = 4;
 
 if (e5) {
   ents[4] = e5;
   count = 5;
 } 
 if (e6) {
   ents[5] = e6;
   count = 6;
 } 

 /* if ents are vertices get the vertices of the face
 */
 int Tface = 2;
 int Tedge = 1;
 int Tvertex =0;

 if(type>Tface)
   return((pRegion)nullptr);

 if( type == Tvertex )
 {
   regions = V_regions((pVertex) ents[0]);
   /* For each face connected to ents[0] see if the other vertices are 
      connected to the face
   */
   
   for (;(region = (pRegion)PList_next(regions, &tmp));) 
   {
     entities = R_vertices(region, 1);
     for( i = 1; i < count; i++)
       if (!PList_inList(entities, ents[i]))
	 break;
     PList_delete(entities);
     if (i == count){
       PList_delete(regions);/* Face is connected to all entities */
       return region;
     }
   }
 }
 else if (type == Tedge)
 {
   regions = E_regions((pEdge) ents[0]);
   /* For each face connected to ents[0] see if the other edges are 
      connected to the face
   */
   
   for (;(region = (pRegion)PList_next(regions, &tmp));) {
     entities = R_edges(region, 1);
     for( i = 1; i < count; i++)
       if (!PList_inList(entities, ents[i]))
	 break;
     PList_delete(entities);
     if (i == count) {
       PList_delete(regions);/* Face is connected to all entities */
       return region;
     }
   }
 }
 
 else if (type == Tface)
 {
   regions = F_regions((pFace) ents[0]);
   /* For each face connected to ents[0] see if the other edges are 
      connected to the face
   */
   
   for (;(region = (pRegion)PList_next(regions, &tmp));) {
     entities = R_faces(region, 1);
     for( i = 1; i < count; i++)
       if (!PList_inList(entities, ents[i]))
	 break;
     PList_delete(entities);
     if (i == count) {
       PList_delete(regions);/* Face is connected to all entities */
       return region;
     }
   }
 } 

 PList_delete(regions);

 return (pFace)nullptr ;
}

