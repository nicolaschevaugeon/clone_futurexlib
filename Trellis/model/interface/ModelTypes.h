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
#ifndef H_ModelTypes
#define H_ModelTypes

//#include "GEntity.h"

/* Geometry type codes */ 
enum GEOM_type {
  GEOM_Unknown,
  GEOM_Point,
  GEOM_Line,
  GEOM_Circle,
  GEOM_Ellipse,
  GEOM_ParametricCurve,
  GEOM_Plane,
  GEOM_Nurb,
  GEOM_Cylinder,
  GEOM_Sphere,
  GEOM_Cone,
  GEOM_Torus,
  GEOM_ParametricSurface
};

/* representation types, returned by GEN_repType */
enum REP_type{
  Unknown,
  Rep_Parametric,
  Rep_Mesh
};

/* code for parametric space status */
enum paramStatus { PAR_UNDEF = -1,
                   PAR_CONT  = 0,
                   PAR_PRBC  = 1, 
                   PAR_PRBNC = 2 
                 };
enum paramType   { PAR_UNKNOWN = -1,
                   PAR_RDEG   = 0,
                   PAR_RNDEG  = 1, 
                   PAR_NRDEG   = 2,
                   PAR_NRNDEG  = 3 
                 };


enum M_ErrorCode { MWARN=3000,
		   MERROR,
		   MFATAL,
		   MMESG
                 };

typedef enum gType{
  Gvertex,
  Gedge,
  Gface,
  Gregion,
  Gvertexuse,
  Gedgeusepair,
  Gloopuse,
  Gfaceuse,
  Gshell
} gType;


#ifdef __cplusplus

template<class T> class SSListCIter; 

extern "C"
{
  typedef class SGModel * pGModel;
//  typedef class AOMD::GEntity * pGEntity;
namespace AOMD{
  typedef class GEntity * pGEntity;
}
  typedef class GRegion * pGRegion;
  typedef class GFace * pGFace;
  typedef class GEdge * pGEdge;
  typedef class GVertex * pGVertex;
  typedef struct PList * pPList;
  
  typedef class GFaceUse * pGFaceUse;
  typedef class GLoopUse * pGLoopUse;
  typedef class GEdgeUsePair * pGEdgeUsePair;
  typedef class GVertexUse * pGVertexUse;
  typedef class GLUIter * pGLUIter;
  typedef class GFUIter * pGFUIter;
  
  typedef SSListCIter<pGRegion> *GRIter;
  typedef SSListCIter<pGFace> *GFIter;
  typedef SSListCIter<pGEdge> *GEIter;
  typedef SSListCIter<pGVertex> *GVIter;
}

#else   /* not __cplusplus */

typedef struct SGModel * pGModel;
typedef struct GEntity * pGEntity;
typedef struct GRegion * pGRegion;
typedef struct GFace * pGFace;
typedef struct GEdge * pGEdge;
typedef struct GVertex * pGVertex;
typedef struct PList * pPList;

typedef struct GFaceUse * pGFaceUse;
typedef struct GLoopUse * pGLoopUse;
typedef struct GEdgeUsePair * pGEdgeUsePair;
typedef struct GVertexUse * pGVertexUse;
typedef struct GLUIter * pGLUIter;
typedef struct GFUIter * pGFUIter;

typedef struct ModelIter *GRIter;
typedef struct ModelIter *GFIter;
typedef struct ModelIter *GEIter;
typedef struct ModelIter *GVIter;

#endif /* not __cplusplus */


#endif

