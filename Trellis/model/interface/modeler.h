/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#ifndef H_modeler
#define H_modeler

#include "ModelTypes.h"
   
#ifdef __cplusplus
extern "C" {
#endif

using AOMD::pGEntity;

/* model operators */

void GM_registerShapes(void);
void GM_registerParasolid(void);

pGModel GM_createFromFile(const char * name);
pGModel GM_createTopoModel(pGModel model);
void GM_delete(pGModel model);

void GM_writeSMD(pGModel model, char *name);

pGEntity GM_entityByTag(pGModel model, int type, int tag);
pGEntity GM_entityByID(pGModel m,int type, int index);
pGEntity GM_entity(pGModel model, int type, int index);

int GM_numVertices(pGModel);
int GM_numEdges(pGModel);
int GM_numFaces(pGModel);
int GM_numRegions(pGModel);

GRIter GM_regionIter(pGModel);
GFIter GM_faceIter(pGModel);
GEIter GM_edgeIter(pGModel);
GVIter GM_vertexIter(pGModel);

pGRegion GRIter_next(GRIter);
pGFace GFIter_next(GFIter);
pGEdge GEIter_next(GEIter);
pGVertex GVIter_next(GVIter);

void GRIter_delete(GRIter);
void GFIter_delete(GFIter);
void GEIter_delete(GEIter);
void GVIter_delete(GVIter);

void GRIter_reset(GRIter);
void GFIter_reset(GFIter);
void GEIter_reset(GEIter);
void GVIter_reset(GVIter);


/*  pGVertex GM_nextVertex(pGModel, void **last); */
/*  pGEdge GM_nextEdge(pGModel, void **last); */
/*  pGFace GM_nextFace(pGModel, void **last); */
/*  pGRegion GM_nextRegion(pGModel, void **last); */
/*  void GM_nextVertexDelete(pGModel, void **last); */
/*  void GM_nextEdgeDelete(pGModel, void **last); */
/*  void GM_nextFaceDelete(pGModel, void **last); */
/*  void GM_nextRegionDelete(pGModel, void **last); */

void GM_bounds(pGModel, double *min, double *max);

pGRegion GM_outerRegion(pGModel);

double GM_tolerance(pGModel model);

/*** model entity operators ***/
int GEN_type(pGEntity);  
int GEN_tag(pGEntity);
int GEN_id(pGEntity);  
int GEN_repType(pGEntity);

void GEN_bounds(pGEntity e, double *min, double *max);
int GEN_inClosure(pGEntity target, pGEntity ent);
int GEN_point(pGEntity ent, double *entPar, double *xyz);
void GEN_reparam(pGEntity target, pGEntity ent, double *entPar, int entDir,
		 double *targetPar);

void GEN_attachDataP(pGEntity , char *tag, void *data);
void * GEN_dataP(pGEntity , char *tag);
int GEN_modifyDataP(pGEntity, char *tag, void * data);
void GEN_attachDataI(pGEntity , char *tag, int data);
int GEN_dataI(pGEntity , char *tag);
int GEN_modifyDataI(pGEntity, char *tag, int data);
void GEN_removeData(pGEntity , char *tag);

void * GEN_getNativePtr(pGEntity);
int GEN_getNativeInt(pGEntity);

/*** regions ops ***/
int GR_inClosure(pGRegion r, pGEntity ent);
pPList GR_faces(pGRegion r);
pPList GR_edges(pGRegion r);
pPList GR_vertices(pGRegion r);

/*** face ops ***/
int GF_inClosure(pGFace f, pGEntity ent);
pPList GF_regions(pGFace f);
pGFaceUse GF_use(pGFace f, int dir);
pPList GF_edges(pGFace f);
pPList GF_vertices(pGFace f);
pGRegion  GF_region(pGFace f, int dir);

void GF_parRange(pGFace f, int axis, double *low, double *high);
int GF_periodic(pGFace f, int dir);

int GF_point(pGFace f, double *par, double *xyz);
void GF_normal(pGFace f, double *par, double *xyz);
void GF_edgeReparam(pGFace f, pGEdge e, double epar, int edir, double *fpar);
void GF_vertexReparam(pGFace f, pGVertex v, double *fpar);

pPList GF_splitU(pGFace f, double u);
pPList GF_splitV(pGFace f, double v);

/*** face use ops ***/
pPList GFU_loopUses(pGFaceUse f);
pGFUIter GFU_loopIter(pGFaceUse f);
int GFUIter_next(pGFUIter fuIter, pGLoopUse *lu);
pGRegion GFU_region(pGFaceUse fu);
pGFace GFU_face(pGFaceUse fu);

/*** loop use ops ***/
pGLUIter GLU_edgeIter(pGLoopUse lu);
int GLUIter_next(pGLUIter luIter, pGEdgeUsePair *edge, int *dir);

/*** edge ops ***/
int GE_inClosure(pGEdge e, pGEntity ent);
pPList GE_regions(pGEdge e);
pPList GE_faces(pGEdge e);
pPList GE_vertices(pGEdge e);
pGVertex GE_vertex(pGEdge e, int dir);
void GE_parRange(pGEdge e, double *low, double *high);
int GE_isSeam(pGEdge e, pGFace face); 
/* note: GE_isSeam will sometimes return true when the edge isn't a seam */
int GE_point(pGEdge e, double par, double *xyz);
double GE_vertexReparam(pGEdge e, pGVertex v);
pPList GE_uses(pGEdge e);

/* depreciated */
int GE_numFaces(pGEdge);
pGFace GE_face(pGEdge,int);

/*** edge use ops ***/
pGEdge GEUP_edge(pGEdgeUsePair eu);
pGVertexUse GEUP_vertexUse(pGEdgeUsePair eu, int dir);
pGFaceUse GEUP_faceUse(pGEdgeUsePair eu, int dir);

/*** vertex ops ***/
pPList GV_regions(pGVertex v);
pPList GV_faces(pGVertex v);
pPList GV_edges(pGVertex v);
void GV_point(pGVertex v, double *xyz);
pPList GV_uses(pGVertex v);

/*** vertex use ops ***/
pGVertex GVU_vertex(pGVertexUse vu);
pPList GVU_edgeUses(pGVertexUse vu);
pPList GVU_faceUses(pGVertexUse vu);

/* depreciated */
pPList GV_numEdges(pGVertex v);
pGEdge GV_edge(pGVertex, int);


/* operators from old interface, all should be emulated in new interface */

/* these operators do not yet have replacements in new interface */
typedef double VECD[3];
int C_deriv(int type, void *tag, double *par, int order, double *deriv);
int C_parStatus(int type, void *tag, int dim);
int C_parType(int type, void *tag, int dim);
int C_parType2(int type, void *tag, int dim,int *nbr, double par[2]);
int C_parOrient(int type, void *tag);
int C_F_edgeDir(void *ftag, void *etag);
int C_tptype(void *tag);  /* Topological entity type */
int C_gmtype(int type, void *tag); /* Geometric entity type */
int C_param(int type, void *tag, double *xyz, double *par);
int C_closestPoint(int type, void *tag, double *inxyz, int seedflag,
                    double *seedxyz, double *seedpar, double *outxyz,
                    double *outparm);
int C_modtol(int type, void *tag, double *tol) ;
int C_gtmten(int type, void * tag1, void ** tag2, double transform[][4]);
int C_endModeler(char *mtype, char *mname) ;
int C_getunv(double *, double *);
int C_lineFaceX(void *ftag, double *lstart, double *lend,
             int max_inter,int *nint, VECD *inter_pts,
             VECD *params, void **boundEntity, int *boundType );
int C_chkContainment(int type, void *ent, int iflag, double *pnt);
int C_chkPtOnBoundary(int type, void *ent, int iflag, double *pnt, int *subtyp,
		      void **subent);
int C_R_useFace(void *region, void *face);


double GEN_tolerance(pGEntity e);
int GF_geomDirection(pGFace e);
int GE_periodic(pGEdge e);
int GF_paramDegeneracies(pGFace f, int dir, double *par);


/* PList  Operators */
pPList  PList_new(void);
pPList  PList_newCopy(pPList );
void PList_delete(pPList );
pPList  PList_copy(pPList ,pPList source);
void PList_clear(pPList);
pPList  PList_appPListUnique(pPList , pPList source);
pPList  PList_appUnique(pPList , void * item);
pPList  PList_append(pPList , void * item);
int PList_size(pPList ); 
void *  PList_item(pPList , int n);
void * PList_next(pPList, void **restart);
int PList_inList(pPList , void * itme);
void PList_remItem(pPList y, void * item);

/* Error Handler */
void M_errHndl(char *msg, char *fct, int severity);

/* Error reporting */
int C_failure(int ifail);

/* Some other functionality used in above routine */
int C_chk_dist( double *point1, double *point2, double tolren );
int C_range_chk( double *entbox1,double *entbox2, int dimension);

/* mappings between old and new operators 
   C_parRange ->  F_parRange, E_parRange 
   C_loadModel -> GM_createFromFile
   C_E_FaceList ->  GE_faces 
   C_E_vertex ->  GE_vertex
   C_F_EdgeList -> GF_edges
   C_V_EdgeList -> GV_edges 
   C_V_FaceList -> GV_faces
   C_onClosure ->  GR_inClosure, GF_inClosure, GE_inClosure 
   C_F_Region ->  GF_region 
   ObjectByTag -> GM_entityByTag 
   TagByObject -> GEN_tag 
   C_point -> GV_point, GE_point, GF_point
   C_normal -> GF_normal
   C_reparam -> GF_edgeReparam, GF_vertexReparam, GE_vertexReparam
   */

#ifdef GM_OLD_OPS
/* these operators have equivalents in the new interface */
int C_reparam(int type1, void *tag1, double *par1,
              int type2, void *tag2, double *par2);
int C_normal(int type, void *tag, double *par, double *xyz);
int C_point(int type, void *tag, double *par, double *xyz);
int C_parRange(int type, void *tag, int dim, double *low, double *high);
int C_loadModel(char *mtype, char *mname);
pPList C_E_FaceList(void *etag);  
void * C_E_vertex(void *etag, int n);
pPList C_F_EdgeList(void *ftag);
pPList C_V_EdgeList(void *vtag);
pPList C_V_FaceList(void *vtag);
int C_onClosure(int type1, void*ent1, int type2, void *ent2);
void *C_F_Region(void *, int ) ; 

/*these operators shouldn't be used anymore*/
void * ObjectByTag(int type, int tag);
int TagByObject(void *obj); 

#endif

#ifdef __cplusplus
}
#endif

#endif
