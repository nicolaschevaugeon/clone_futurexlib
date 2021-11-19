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
#ifndef H_AOMD
#define H_AOMD

#ifndef SIM
#include "AOMDfwd.h"
#include "modeler.h"

#ifdef __cplusplus
extern "C" {
  typedef char MeshDataId;

  typedef class AOMD::mMesh  * pMesh;  
  
  typedef class AOMD::mEntity * pRegion;
  typedef class AOMD::mEntity * pFace;
  typedef class AOMD::mEntity * pEdge;
  typedef class AOMD::mEntity * pVertex;
  typedef class AOMD::mEntity * pPoint;
  typedef class MNode * pNode;
  
  typedef class MeshList * pMeshList;
  typedef class MeshListIter * pMeshListIter;
  
  typedef class AOMD::mEntity * pEntity;
  
  typedef class AOMD::mIteratorNew *RIter;
  typedef class AOMD::mIteratorNew *FIter;
  typedef class AOMD::mIteratorNew *EIter;
  typedef class AOMD::mIteratorNew *VIter;

//  typedef class AOMD::AOMD_LoadBalancerCallbacks MigrationCB;

#else   /* not __cplusplus */

  typedef struct SMesh * pMesh;
  
  typedef struct Region * pRegion;
  typedef struct Face * pFace;
  typedef struct Edge * pEdge;
  typedef struct Vertex * pVertex;
  typedef struct Point * pPoint;
  typedef struct Node * pNode;

  typedef struct MeshList * pMeshList;
  typedef struct MeshListIter * pMeshListIter;

  typedef void * pEntity;

  typedef struct MeshIter *RIter;
  typedef struct MeshIter *FIter;
  typedef struct MeshIter *EIter;
  typedef struct MeshIter *VIter;

#endif /* not __cplusplus */

  typedef enum eType{
    Tvertex,
    Tedge,
    Tface,
    Tregion
  } eType;


  /*** Initialization, logging, execution control ***/
  void MS_logOn(const char *name);
  void MS_logOff();

  void MS_init(void);
  void MS_exit(void);

  void MS_abort();

  /*** Error handling routines ***/

  typedef enum MeshingErrorFix {
    None,
    Refine,
    ModelProblem
  } MeshingErrorFix;

  typedef enum MeshingErrorType {
    ModelEnt,
    XYZ
  } MeshingErrorType;

#ifdef __cplusplus
  typedef class MeshingError * pMeshingError;
#else
  typedef struct MeshingError * pMeshingError;
#endif  /* __cplusplus */

  pPList MS_getErrorList(pMesh mesh);
  MeshingErrorType MS_errorType(pMeshingError merror);
  MeshingErrorFix MS_errorFix(pMeshingError merror);
  pGEntity MS_errorEnt(pMeshingError merror);
  double * MS_errorXYZ(pMeshingError merror, double *xyz);

  /*** mesh object creation ***/
  pMesh MS_newMesh(pGModel);

  /*** mesh generation ***/
  int MS_generateMesh(pMesh, int unused);
  int MS_generateSurfaceMesh(pMesh, int unused);
  int MS_generateVolumeMesh(pMesh, int unused);

  int MS_curveMesh(pMesh mesh);

  int MS_meshModelFace(pMesh mesh, pGFace face);

  /*** model entity mesh control operators ***/
  void MS_setLocalMeshSize(pMesh, pGEntity, int type, double size);
  void MS_setLocalMeshSizeExp(pMesh mesh, pGEntity ent, int type, const char *sizeexp);
  void MS_setGlobalMeshSize(pMesh, int type, double size);
  void MS_setLocalMinCurvSize(pMesh mesh, pGEntity ent, int type, double size);
  void MS_setGlobalMinCurvSize(pMesh mesh, int type, double size);
  void MS_setLocalMeshCurv(pMesh, pGEntity ent, int type, double par);
  void MS_setGlobalMeshCurv(pMesh, int type, double par);

  void MS_setMeshMatch(pMesh, pGEntity src, pGEntity dest, double *trans,
		       double *axis, double *rpnt, double theta);
  void MS_setEdgeMesh(pMesh mesh, pGEdge gedge, int np, double *pars, double *xyzs);
  void MS_setBoundaryLayer(pMesh mesh, pGEntity ent, int side,
			   int type, double t0, double T, int nLayer);

  void MS_setNoMesh(pMesh, pGEntity);

  /*** mesh generation options */
  void MS_setSurfaceMeshType(pMesh, int type);
  void MS_setBLMeshType(pMesh mesh, int type);
  void MS_setVolumeMeshType(pMesh, int type);

  void MS_setSnapWhenMatching(pMesh, int onoff);

  void MS_setMaxEntities(pMesh, int maxv, int maxe, int maxf, int maxr);

  /********************/
  /* Mesh Operators  */
  /********************/
  void M_delete(pMesh);
  int M_load_old(pMesh, const char *filename);
  int M_load(pMesh, const char *filename);
  pGModel M_model(pMesh mesh);

  int M_numRegions(pMesh);
  int M_numFaces(pMesh);
  int M_numEdges(pMesh);
  int M_numVertices(pMesh);
  int M_numPoints(pMesh);

  pRegion M_region(pMesh, int n);
  pFace M_face(pMesh, int n);
  pEdge M_edge(pMesh, int n);
  pVertex M_vertex(pMesh, int n);

  RIter M_regionIter(pMesh mesh);
  FIter M_faceIter(pMesh mesh);
  EIter M_edgeIter(pMesh mesh);
  VIter M_vertexIter(pMesh mesh);

  void M_writeSMS(pMesh mesh, const char *name, int version);

  RIter M_classifiedRegionIter(pMesh mesh, pGEntity ent);
  FIter M_classifiedFaceIter(pMesh mesh, pGEntity ent, int closure);
  EIter M_classifiedEdgeIter(pMesh mesh, pGEntity ent, int closure);
  VIter M_classifiedVertexIter(pMesh mesh, pGEntity ent, int closure);
  pVertex M_classifiedVertex(pMesh mesh, pGVertex ent);

  int M_numClassifiedRegions(pMesh mesh, pGEntity ent);
  int M_numClassifiedFaces(pMesh mesh, pGEntity ent);
  int M_numClassifiedEdges(pMesh mesh, pGEntity ent);
  int M_numClassifiedVertices(pMesh mesh, pGEntity ent);

  /*************************/
  /* Refinement  Operators */
  /*************************/
  void MS_setPointRefinement(pMesh mesh, double size, double pt[3]);
  void MS_setLineRefinement(pMesh mesh, double size, double start[3], double end[3]);
  void MS_setPolyLineRefinement(pMesh mesh, double size, int numpts, double **pts);
  void MS_setPlaneRefinement(pMesh mesh, double size, double center[3], double wdir[3], double hdir[3]);
  void MS_setCircleRefinement(pMesh mesh, double size, double radius, double center[3], double normal[3]);
  void MS_setCylinderRefinement(pMesh mesh, double size, double radius, double length, double center[3], double normal[3]);
  void MS_setCubeRefinement(pMesh mesh, double size, double center[3], double wdir[3], double hdir[3], double ddir[3]);
  void MS_setSphereRefinement(pMesh mesh, double size, double radius, double center[3]);

  /********************/
  /* Entity Iter Ops  */
  /********************/

  pRegion RIter_next(RIter);
  void RIter_delete(RIter);
  void RIter_reset(RIter);
  pFace FIter_next(FIter);
  void FIter_delete(FIter);
  void FIter_reset(FIter);
  pEdge EIter_next(EIter);
  void EIter_delete(EIter);
  void EIter_reset(EIter);
  pVertex VIter_next(VIter);
  void VIter_delete(VIter);
  void VIter_reset(VIter);



  /********************/
  /* Entity Operators */
  /********************/
  void EN_setID(pEntity , int id);
  int EN_id(pEntity );

  int EN_type(pEntity);
  int EN_whatInType(pEntity);
  pGEntity EN_whatIn(pEntity);

  pPList EN_getMatchingEnts(pEntity ent, pGEntity filter);

  /********************/
  /* Region Operators */
  /********************/

  int R_numFaces(pRegion );  /* gets number of faces on a region */
  pFace R_face(pRegion , int n);  /* gets n'th face on region */
  int R_faceDir(pRegion , int n);
  pGRegion R_whatIn(pRegion );  /* gets classification of region */
  int R_whatInType(pRegion);
  pPList R_vertices(pRegion, int);
  pPList R_edges(pRegion, int );
  pPList R_faces(pRegion, int );
  int R_inClosure(pRegion, pEntity ent);
  int R_dirUsingFace(pRegion , pFace face);

  /********************/
  /* Face Operators */
  /********************/

  pGEntity  F_whatIn(pFace );  /* gets classification of face */
  int F_whatInType(pFace);
  int F_numEdges(pFace  );     /* get number of edges around face */
  pEdge F_edge(pFace  , int n); /* get edge n on face */
  int F_edgeDir(pFace  , int n);
  pPList F_edges(pFace , int dir,pVertex vtx);
  pPList F_vertices(pFace, int dir);
  pVertex F_vertex(pFace pf, int n);
  pRegion F_region(pFace, int );
  int F_inClosure(pFace, pEntity ent);
  int F_dirUsingEdge(pFace, pEdge edge);
  pPList F_regions(pFace);
  int F_numRegions(pFace);

  /********************/
  /* Edge Operators */
  /********************/
  pGEntity  E_whatIn(pEdge );  /* gets classification of an edge*/
  int E_whatInType(pEdge);
  int E_numFaces(pEdge  );     /* get the number of faces using an edge */
  int E_numRegions(pEdge);
  pFace E_face(pEdge, int n);
  pVertex E_vertex(pEdge  , int n); /* get vertex n on face */
  pPList E_vertices(pEdge); /* get vertex n on face */
  
  int E_numPoints(pEdge);
  pPoint E_point(pEdge , int n);
  pFace E_otherFace(pEdge , pFace, pRegion);
  pVertex  E_otherVertex(pEdge , pVertex v);
  int E_inClosure(pEdge, pEntity ent);
  pPList E_regions(pEdge);

  /********************/
  /* Point Operators  */
  /********************/
  double P_x(pPoint );
  double P_y(pPoint );
  double P_z(pPoint );

  void P_setID(pPoint , int id);
  int P_id(pPoint);

  double P_param1(pPoint);
  void P_param2(pPoint, double *p1, double *p2, int *ptch);

  /********************/
  /* Vertex Operators */
  /********************/

  pGEntity  V_whatIn(pVertex );  /* gets classification of vertex */
  int V_whatInType(pVertex);
  int V_numEdges(pVertex  );     /* gets number of edges using vertex */
  int V_numFaces(pVertex);
  int V_numRegions(pVertex);
  pEdge V_edge(pVertex  , int n); /* get edge n on vertex */
  pPList V_edges(pVertex); /* get edge n on vertex */
  pPList V_regions(pVertex);
  pPList V_faces(pVertex );
  pPoint  V_point(pVertex);
  void V_coord(pVertex, double *xyz);
  void V_setPoint(pVertex, pPoint);


  /*** Message Handling ***/
  typedef enum MS_MsgLevel {
    MS_InfoMsg=0,
    MS_DebugMsg,
    MS_WarningMsg,
    MS_ErrorMsg,
    MS_InternalErrorMsg
  } MS_MsgLevel;

  typedef void (MS_MESSAGEHANDLER)(MS_MsgLevel type,const char *file, int line, const char *msg);

  void MS_setMessageHandler(MS_MESSAGEHANDLER *);

  /*** Progress Handling ***/
  typedef void (MS_PROGRESSHANDLER)(const char *what, int level, int startVal, int endVal, int currentVal, void *);

  void MS_setProgressHandler(MS_PROGRESSHANDLER *);

  typedef unsigned int MarkID;
  unsigned int EN_mark(pEntity entity, MarkID id);
  void EN_setMark(pEntity entity, MarkID id, unsigned int value);

pPList E_faces(pEdge e);
pEdge E_exist(pVertex v1, pVertex v2);
pFace F_exist(int type, pEntity e1, pEntity e2, pEntity e3, pEntity e4);
pRegion R_exist(int type, pEntity e1, pEntity e2, pEntity e3, 
                         pEntity e4, pEntity e5, pEntity e6);

int M_findE(pMesh mesh, pEntity ent, int i);

#ifdef __cplusplus
}
#endif

#endif	/* ifndef SIM */
#endif /* ifndef H_AOMD */


