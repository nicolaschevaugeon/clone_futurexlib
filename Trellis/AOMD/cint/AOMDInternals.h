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
#ifndef H_AOMDInternal
#define H_AOMDInternal

#include "AOMD.h"
#ifndef SIM
#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
  // SMALL DIFFERENCE WITH MESHSIM !!!!
typedef unsigned int pMeshDataId;
typedef class AttachableData *pAttachableData;
typedef class AttachDataId *pAttachDataId;

#else  /* not _cplusplus */

typedef struct MeshDataId * pMeshDataId;
typedef struct AttachDataId *pAttachDataId;
typedef struct AttachableData *pAttachableData;

#endif /* not _cplusplus */

/* these must match the callback types in MeshDataId.h */
enum { CBdelete, CBmigrateOut, CBmigrateIn };
typedef int (*CBfunc)(pAttachableData, pAttachDataId, int, void**, void*);

  /* these data operators are obsolete and retained for compatibility */
  /* new code should use the new operators defined below */
void EN_attachDataP(pEntity , const char *tag, void *data);
void * EN_dataP(pEntity , const char *tag);
int EN_modifyDataP(pEntity, const char *tag, void * data);

void EN_attachDataI(pEntity , const char *tag, int data);
int EN_dataI(pEntity , const char *tag);
int EN_modifyDataI(pEntity, const char *tag, int data);

void EN_removeData(pEntity , const char *tag);

  /* new data operators */
pMeshDataId MD_newMeshDataId(const char *tag);
pMeshDataId MD_lookupMeshDataId(const char *tag);
void MD_deleteMeshDataId(pMeshDataId id);
void MD_setMeshCallback(pMeshDataId id, int event, CBfunc f, void *cbdata);
void MD_removeMeshCallback(pMeshDataId id, int event);

void EN_attachDataInt(pEntity ent, pMeshDataId id, int value);
void EN_attachDataDbl(pEntity ent, pMeshDataId id, double value);
void EN_attachDataPtr(pEntity ent, pMeshDataId id, void * value);

void EN_deleteData(pEntity ent, pMeshDataId id);

int EN_getDataInt(pEntity ent, pMeshDataId id, int *value);
int EN_getDataDbl(pEntity ent, pMeshDataId id, double *value);
int EN_getDataPtr(pEntity ent, pMeshDataId id, void **value);

void EN_modifyDataInt(pEntity ent, pMeshDataId id, int value);
void EN_modifyDataDbl(pEntity ent, pMeshDataId id, double value);
void EN_modifyDataPtr(pEntity ent, pMeshDataId id, void * value);

int EN_inClosure(pEntity, pEntity);

  /* mesh entity creation routines */
pRegion M_createR(pMesh mesh, int nFace, pFace *faces, int *dirs, pGEntity gent);
pFace M_createF(pMesh mesh, int nEdge, pEdge *edges, int *dirs, pGEntity gent);
pEdge M_createE(pMesh mesh, pVertex v1, pVertex v2, pGEntity ent);
pVertex M_createVP(pMesh mesh, double x, double y, double z, double *param,
		    int patch, pGEntity ent);
pVertex M_createVP2(pMesh mesh, double *xyz, double *param,
                     int patch, pGEntity ent);
pPoint M_createP(pMesh mesh, double x, double y, double z, double *param,
		   int unused, pGEntity ent);

  /* mesh entity deletion routines */
void M_removeRegion(pMesh, pRegion region);
void M_removeFace(pMesh, pFace face);
void M_removeEdge(pMesh, pEdge edge);
void M_removeVertex(pMesh, pVertex vertex);
void M_removePoint(pMesh, pPoint point);

  /* point routines */
pPoint P_new(void);
void P_delete(pPoint);

void P_setPos(pPoint  , double x, double y, double z);
  /// AOMD SPECIFIC 
void AOMD_P_setParametricPos(pPoint  , double x, double y, double z);
void P_setParam1(pPoint, double param);
void P_setParam2(pPoint, double p1, double p2, int ptch);

  /* obsolete iterator functions */
pRegion M_nextRegion(pMesh, void **restart);
pFace M_nextFace(pMesh, void **restart);
pEdge M_nextEdge(pMesh, void **restart);
pVertex M_nextVertex(pMesh, void **restart);
pPoint M_nextPoint(pMesh, void **restart);

pRegion M_nextRegionCancel(pMesh, void **restart);
pFace M_nextFaceCancel(pMesh, void **restart);
pEdge M_nextEdgeCancel(pMesh, void **restart);
pVertex M_nextVertexCancel(pMesh, void **restart);
pPoint M_nextPointCancel(pMesh, void **restart);

int middlePoint
(pEdge,     // an edge to be split
 double*,   // the boundary location of the vertex created in splitting
 double*);  // parameter of the vertex in case parametric space exists



/* Change with meshsim, one extra parameter */

void V_merge2(pMesh pM, pVertex, pVertex);
void V_merge(pVertex, pVertex);
void E_merge2(pMesh pM, pEdge, pEdge);
void E_merge(pEdge, pEdge);
void F_merge2(pMesh pM, pFace, pFace);
void F_merge(pFace, pFace);

void EN_setWhatIn (pMesh pm, pRegion e, pGEntity what);
void R_setWhatIn(pRegion  , pGEntity what);
void F_setWhatIn(pFace  , pGEntity what);
void F_chDir(pFace);
void E_setWhatIn(pEdge  , pGEntity what);
void V_setWhatIn(pVertex, pGEntity what); /*sets classification of vertex */

void E_setPoint(pEdge , pPoint pt);
void V_setPoint(pVertex, pPoint pt);

void PList_app2Lst(pPList*, pPList, int);

pPList V_bdryFaces(pVertex);
int F_checkFlat(pFace, double*);

void V_setSize(pVertex,double);
double V_size (pVertex pe);

  /** MARKS (should NOT be used)*/
  
  MarkID MD_registerMark(char *name, unsigned int size, unsigned int flags);
  MarkID MD_markID(char *name);
  unsigned int EN_mark(pEntity entity, MarkID id);
  void EN_setMark(pEntity entity, MarkID id, unsigned int value);


  /*
    return sense of edge use in face
    0 = error, model not loaded or edge does not bound face
    1 = same,
    -1 = opposite
    2 = used in both directions
    -2 = unknown
    ftag = a pointer to model face
    etag = a pointer to model edge
  */
  int C_F_edgeDir(void *ftag, void *etag);


#ifdef __cplusplus
}
#endif
#endif
#endif /* not H_AOMDInternal */
