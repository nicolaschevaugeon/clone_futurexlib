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
#include <cstdlib>
#include <cstring>
#include "SGModel.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
#include "GLoopUse.h"
#include "GShell.h"
#include "GRegion.h"
#include "GVPoint.h"
#include "GEPoint.h"
#include "GFPoint.h"
#include "GEdgeUsePair.h"
#include "SPoint3.h"
#include "Pair.h"
#include "Range.h"
#include "SSList.h"
#include "toPList.h"
#include "modeler.h"

#include <cfloat>
#define MAXDOUBLE DBL_MAX

using std::cerr;
//
// C- callable operators to methods on model object
//
// Right now can deal with a single model at a time
//

SGModel *theModel=nullptr;       // refers to a generic geo. modeler instance


#ifdef __cplusplus
extern "C"
{
#endif


//
// Load the geometric model
// 
// mtype = modeler name : "shapes","parasolid"
// mname = model file name : "model.geom","model"
//
int C_loadModel(char *, char *)
{
  
  if(!theModel)
    cerr << "C_loadModel - Model is not loaded\n";
  //theModel = SGModel::createModel(mname,mtype);

  return 1;
}

//
// return the modeler object pointer represented by this tag
// type : 0=vertex, 1=edge, 2=face, 3=region
//
void * ObjectByTag(int type, int tag)
{
  AOMD::GEntity *obj ;

  switch( type )
  {
    case 0 :  // vertex
     obj = theModel->vertexByTag(tag);
     break;
    case 1 :  // edge
     obj = theModel->edgeByTag(tag);
     break;
    case 2 :  // face
     obj = theModel->faceByTag(tag);
     break;
    case 3 :  // region
     obj = theModel->regionByTag(tag);
     break;
    default :
     obj = nullptr ;
     break ;
  }
  return (void *)obj ;
}

//
// return unique modeler tag for this modeler object 
//
int TagByObject(void *obj)
{
  if( obj == (void *)nullptr )
   return 0 ;

  AOMD::GEntity *gobj = (AOMD::GEntity *)obj ;
  return gobj->tag() ;
}

//
// Return xyz on a model entity at a given parameter location
//
// type = type of model entity (TopoType::Vertex, TopoType::Edge, TopoType::Face)
// tag  = model entity tag ( got using EN_tag() )
// par  = parameter value (1 for TopoType::Edge, 2 for TopoType::Face)
// xyz  = returned postion vector
//
int C_point(int type, void *tag, double *par, double *xyz)
{
  AOMD::SPoint3 pnt ;
  int valid = 1;

  switch(type)
  {
    case 0 : // vertex
    {
      GVertex *vtx = (GVertex *)tag ;
      pnt = vtx->point();
      break ;
    }
    case 1 : // edge
    {
      GEdge *edg = (GEdge *)tag ;
      GEPoint epnt = edg->point(par[0]);
      valid = epnt.isValid();
      pnt = epnt;
      break ;
    }
    case 2 : // face
    {
      GFace *fac = (GFace *)tag ;
      GFPoint fpnt = fac->point(par[0],par[1]);
      valid = fpnt.isValid();
      pnt = fpnt;
      break ;
    }
    default : // error
      return 0;
  }
  xyz[0] = pnt.x() ; 
  xyz[1] = pnt.y(); 
  xyz[2] = pnt.z() ;
  return valid;
}

//
// reparameterize a point wrt to higher order model entity
// type1 = type of lower order model entity
// tag1  = pointer to lower order model entity
// type2 = type of higher order model entity
// tag2  = pointer to higher order model entity
// par1  = parameters wrt tag1
// par2  = parameters wrt tag2 (output)
//
int C_reparam(int type1, void *tag1, double *par1,
              int type2, void *tag2, double *par2)
{
  // given par1 wrt tag1 of type1, get par2 wrt tag2 of type2

  // if type1 is higher order than type2 then return w/ error
  if( type1 > type2 ) 
    return -1 ;

  // if the types are of the same order then the tags must match
  if( type1 == type2 )
  {
    if( TagByObject(tag1) != TagByObject(tag2) )
      return -1 ;
    else
    {
      par2[0] = par1[0]; 
      par2[1] = par1[1]; 
      par2[2] = par1[2]; 
      return 1 ;
    }
  }

  SPoint2 pt ;
  GVPoint vPoint ;
  GEPoint ePoint ;

  if( type1 == TopoType::Vertex ){
    GVertex *inVtx = (GVertex *)tag1 ;
    vPoint = inVtx->point();
  }else if( type1 == TopoType::Edge ){
    GEdge *inEdg = (GEdge *)tag1 ;
    ePoint = inEdg->point(par1[0]);
  }
 
  if( type2 == TopoType::Edge ){
    GEdge *outEdg = (GEdge *)tag2 ;
    if( type1 == TopoType::Vertex )
      par2[0] = outEdg->param(vPoint);
  } else if( type2 == TopoType::Face ) {
    GFace *outFac = (GFace *)tag2 ;
    if( type1 == TopoType::Vertex ){
      pt = outFac->param(vPoint);
      par2[0] = pt[0]; par2[1] = pt[1];
    }else if (type1 == TopoType::Edge){
      pt = outFac->param(ePoint);
      par2[0] = pt[0]; par2[1] = pt[1];
    }
    par2[2] = (double)outFac->tag();
  }
  return 1;
}

//
// return modeler derivatives wrt a given model entity
// type = type of model entity
// tag  = model entity pointer
// par  = parameter at which to evaluate derivs.
// order= derivative order
// dvec = derivative components returned as a vector. in case of
//        surfaces contains derivatives wrt 'u' followed by deriv.
//        wrt 'v'.
//
int C_deriv(int type, void *tag, double *par, int order, double *dvec)
{
  int retval = 1 ;
  if( type == TopoType::Edge )
  {
    if(order == 1 ){
      GEdge *edg = (GEdge *)tag ;
      SVector3 vec1 = edg->firstDer(par[0]);
      dvec[0] = vec1.x() ;
      dvec[1] = vec1.y() ;
      dvec[2] = vec1.z() ;
    } else if (order ==2) {
      GEdge *edg = (GEdge*)tag;
      edg->nthDerivative(par[0],2,dvec);
    } else
      return 0;
  }
  else if (type == TopoType::Face)
  {
    SPoint2 uv(par[0],par[1]);
    GFace *fac = (GFace *)tag ;

    if( order == 1 )
    {
      Pair <SVector3,SVector3> der1 = fac->firstDer(uv);
      SVector3 vec1 = der1.first() ;
      dvec[0] = vec1.x() ;
      dvec[1] = vec1.y() ;
      dvec[2] = vec1.z() ;
      SVector3 vec2 = der1.second() ;
      dvec[3] = vec2.x() ;
      dvec[4] = vec2.y() ;
      dvec[5] = vec2.z() ;
    }
    else if( order == 2 || order == 3 ) {
     
      /* WARNING: using undocumented function within the modeler for
         shapes
      */
      //double dtemp[3],dtemp1[3] ;
      if( !fac->nthDerivative(uv,order,dvec) )
        retval = 0;
    }
    else
      return 0 ;
  }
  return 1 ;
}

//
// Normal on a model face at given parameters
//
int C_normal(int type, void *tag, double *par, double *xyz)
{

  // if model is not loaded, return w/ error
  // also, normal defined only for TopoType::Face
  if (type != TopoType::Face) 
    return -1 ;


  SPoint2 pt(par[0],par[1]);

  GFace *fac = (GFace *)tag ;
  SVector3 norm = fac->normal(pt);

  xyz[0] = norm.x();
  xyz[1] = norm.y();
  xyz[2] = norm.z();

  return 1 ;
}

//
// return type of underlying model geometry
//
int C_gmtype(int type, void *tag)
{
  int gtype ;

  // get the tag-th entity 
  if( type == TopoType::Vertex ){
    GVertex *inVtx = (GVertex *)tag ;
    gtype = inVtx->geomType(); 
  }else if( type == TopoType::Edge ){
    GEdge *inEdg = (GEdge *)tag ;
    gtype = inEdg->geomType(); 
  }else if( type == TopoType::Face ){
    GFace *inFac = (GFace *)tag ;
    gtype = inFac->geomType(); 
  }
  switch(gtype)
  {
    case GeomType::Point :
      return GEOM_Point ; 
    case GeomType::Line  :
      return GEOM_Line ; 
    case GeomType::Circle :
      return GEOM_Circle ; 
    case GeomType::Ellipse :
      return GEOM_Ellipse ; 
    case GeomType::ParametricCurve :
      return GEOM_ParametricCurve ; 
    case GeomType::Plane :
      return GEOM_Plane ; 
    case GeomType::Nurb :
      return GEOM_Nurb ; 
    case GeomType::Cylinder :
      return GEOM_Cylinder ; 
    case GeomType::Sphere :
      return GEOM_Sphere ; 
    case GeomType::Cone :
      return GEOM_Cone ; 
    case GeomType::Torus :
      return GEOM_Torus ; 
    case GeomType::ParametricSurface :
      return GEOM_ParametricSurface ; 
    default :
      return GEOM_Unknown ; 
  }
}

//
// return periodicity information about geometry underlying a model
// entity in a given parameter direction
// dim = 0 (u direction)
//     = 1 (v direction)
//
int C_parStatus(int type, void *tag, int dim)
{
  int pflag{0}, cflag{0} ;

  // get the tag-th entity 
  if( type == TopoType::Edge ){
    GEdge *inEdg = (GEdge *)tag ;
    pflag = inEdg->periodic(dim); 
    cflag = inEdg->continuous(dim); 
  }else if( type == TopoType::Face ){
    GFace *inFac = (GFace *)tag ;
    pflag = inFac->periodic(dim); 
    cflag = inFac->continuous(dim); 
  }
  if( pflag )
  {
    if( cflag )
     return PAR_PRBC ;
    else
     return PAR_PRBNC ;
  }
  else
  {
    if( cflag )
     return PAR_CONT ;
  }
  return PAR_UNDEF ;
}

//
// return type information about geometryparameterization underlying 
// a model entity
// dim = 0 (u direction)
//     = 1 (v direction)
//
int C_parType(int type, void *tag, int dim)
{
  int dflag{0} ;

  // get the tag-th entity 
  if( type == TopoType::Edge ){
    GEdge *inEdg = (GEdge *)tag ;
    dflag = inEdg->degenerate(dim); 
  }else if( type == TopoType::Face ){
    GFace *inFac = (GFace *)tag ;
    dflag = inFac->degenerate(dim); 
  }
  if( dflag ) 
    return PAR_RDEG ;     // degenerate
  else
    return PAR_RNDEG ;    // non-degerate
}

//
// return parameter range of underlying geometry
// dim = 0 (u direction)
//     = 1 (v direction)
// bpar= range start
// epar= range end
//
int C_parRange(int type, void *tag, int dim, double *bpar, double *epar)
{
  if (type < TopoType::Edge || type > TopoType::Face )
    return -1;

  AOMD::Range<double> range ;

  if( type == TopoType::Edge )
  {
    if( dim != 0 )
      return -1 ;
    GEdge *edge = (GEdge *)tag ;
    range = edge->parBounds(dim);
  }
  else if (type == TopoType::Face)
  {
    if( dim != 0 && dim != 1 )
      return -1 ;
    GFace *face = (GFace *)tag ;
    range = face->parBounds(dim);
  }
  *bpar = range.low();
  *epar = range.high();
  return 1 ;
}

//
// return parameteric orientation of the model face wrt the 
// surface parameterisation
// 1 = same orientation
// -1 = opp. orientation
//
int C_parOrient(int type, void *tag)
{
  if( type == TopoType::Edge )
  {
    GEdge *edge = (GEdge *)tag ;
    return edge->geomDirection();
  }
  else if( type == TopoType::Face )
  {
    GFace *face = (GFace *)tag ;
    return face->geomDirection();
  }
  cerr << "Warning: C_parOrient called for other than face or edge\n";
  return 1;
}

//
// return the list of model faces sharing a model edge
// tag = ptr to model edge
//
pPList C_E_FaceList(void *tag)
{
  GEdge *edge = (GEdge *)tag ;

  // get the SSList of faces sharing this edge
  SSList<GFace*> faces(edge->faces());
  SSListIter<GFace*> fiter(faces);
  GFace *face ;
  
  pPList ftags = PList_new();
  while( fiter(face) )
  {     
    PList_append(ftags,(pGEntity)face); 
  }

  return ftags ;
}

//
// get the i-th end vertex of a model edge
// etag = model edge pointer
// n=0 start vertex, 1 = end vertex
//
void * C_E_vertex(void *etag, int n)
{

  GEdge *edge = (GEdge *)etag ;
  GVertex *vtx = edge->vertex(n);
    
  return (void *)vtx;
}

//
// return sense of edge use in face 
// 0 = error, model not loaded or edge does not bound face
// 1 = same, 
// -1 = opposite
// 2 = used in both directions
// -2 = unknown
// ftag = face pointer
// etag = edge pointer
//
int C_F_edgeDir(void *ftag, void *etag)
{
  int sense = 0;
  int nuse = 0;
  GFace *face = (GFace *)ftag ;
  GEdge *inedg = (GEdge *)etag ;

  // locate the edge in the list of edges by getting the loops
  // for this face and then looking at the edges of the loops
  
  SSListCIter<GLoopUse*> LopIter = face->use(1)->firstLoopUse();
  GLoopUse *loop;
  while( LopIter(loop) ){
    SSListCIter<GEdgeUse*> eui = loop->firstEdgeUse();
    GEdgeUse *eu;
    while( eui(eu) ){
      if( eu->edge() == inedg ){
        switch( eu->dir() ){
	case 0 :  // opposite direction
	  sense = ( nuse ? 2 : -1 );
	  nuse++;
	  break ;
	case 1 :  // same direction
	  sense = ( nuse ? 2 : 1 );
	  nuse++;
	  break ;
	default:  // can't happen
	  cerr << "C_F_edgeDir - error\n";
	  break ;
        }
      }
    }
    if( nuse )
      break ;
  }
  
  return sense ;
}

//
// return the list of ptrs of edges sharing a vertex
//
pPList C_V_EdgeList(void *tag)
{
  GVertex *vert = (GVertex *)tag ;
  return toPList(vert->edges());
}

//
// return list of edges bounding a face
//
pPList C_F_EdgeList(void *ftag)
{
  GFace *face = (GFace *)ftag ;
  return toPList(face->edges());
}

pPList C_V_FaceList(void *vtag)
{
  GVertex *vertex = (GVertex*)vtag;
  return toPList(vertex->faces());
}

//
// check if ent1 belongs to closure of ent2 ( return 1 if true, 0 else )
//
int C_onClosure(int type1, void *ent1, int type2, void *ent2)
{
  if( type1 > type2 )
    return 0 ;
  else if( type1 == type2 )
    if( ((AOMD::GEntity *)ent1)->tag() == ((AOMD::GEntity *)ent2)->tag() )
      return 1 ;
    else
      return 0 ;
  else
  {
    if( type2 == TopoType::Edge )
    {
       /* type1 must be a model vertex */
       if( ((AOMD::GEntity *)ent1)->tag() == (((GEdge *)ent2)->vertex(0))->tag() ||
           ((AOMD::GEntity *)ent1)->tag() == (((GEdge *)ent2)->vertex(1))->tag() )
         return 1 ;
       else
         return 0 ;
    }
    else if ( type2 == TopoType::Face )
    {
       /* type1 must be edge for vertex */
        SSListCIter<GLoopUse*> LopIter = ((GFace*)ent2)->use(1)->firstLoopUse();
        GLoopUse *loop ;
        while( LopIter(loop) )
        {
          SSList<GEdge*> edges = loop->edges();
          SSListIter<GEdge*> EdgIter(edges);
          GEdge *edge ;
          while( EdgIter(edge) )
	  {
            if( type1 == TopoType::Edge )
	    {
	      if( edge->tag() == ((AOMD::GEntity *)ent1)->tag() )
                return 1;
            }
            else
	    {
	      if( ((AOMD::GEntity *)ent1)->tag() == (edge->vertex(0))->tag() ||
		  ((AOMD::GEntity *)ent1)->tag() == (edge->vertex(1))->tag() )
                return 1 ;
            }
          }
        }
        return 0;        
    }
    else if ( type2 == TopoType::Region )
     return 0;
  }
  return 0;
}


//int C_param(int type, void *tag, double *xyz, double *par)
//{
//  par[0] = 0.0 ;
//  if( type == TopoType::Region )
//    return 0;
//
//  return 0 ;
//}

/* ---------------------------------------------------------------
   closest point on a model entity for an arbitrary point in space
   not necessarily on the model entity
   
   Input :

     type : TopoType::Edge/TopoType::Face
     tag  : entity object
     inxyz: input point
 seedflag : 0  => NO seed point or parameter
            1  => seed point only
            2  => seed parameter only
            3  => seed point and parameter
 seedxyz  : seed point position
 seedpar  : seed point parameter

   Output :

 outxyz   : closest point
 outparm  : parameters of the closest point
---------------------------------------------------------------- */
int C_closestPoint(int type, void *tag, double *inxyz, int seedflag,
                    double *seedxyz, double *seedpar, double *outxyz,
                    double *outparm)
{
  // create an Spoint3 object
  AOMD::SPoint3 point(inxyz[0],inxyz[1],inxyz[2]);

  if( type == TopoType::Vertex || type == TopoType::Region )
    return 0 ;
  else if( type == TopoType::Edge )
  {
    GEPoint opoint = ((GEdge *)tag)->closestPoint(point);
    outparm[0] = opoint.par() ;
    outxyz[0] = opoint.x() ;
    outxyz[1] = opoint.y() ;
    outxyz[2] = opoint.z() ;
  }
  else if( type == TopoType::Face )
  {
    GFPoint opoint = ((GFace *)tag)->closestPoint(point); 
    SPoint2 par = opoint.par() ;
    outparm[0] = par.x() ;
    outparm[1] = par.y() ;
    outparm[2] = (double) ((AOMD::GEntity *)tag)->tag();
    outxyz[0] = opoint.x() ;
    outxyz[1] = opoint.y() ;
    outxyz[2] = opoint.z() ;
  }

  return 1 ;
}

//
// return sense of face use in region
// -1 = face normal points into region
// 1 = face normal points outof region
// 0 = used both ways, face embedded in region
//
int C_R_faceDir(void *rtag, void *ftag)
{
  GFace *face = (GFace *)ftag ;
  GRegion *region = (GRegion *)rtag ;

  GRegion *r0 = face->region(0);
  GRegion *r1 = face->region(1);

  int sense=2;
  if(r0==region)
    sense = 1;
  if(r1==region)
    sense = -1;
  if(r0==region && r1==region)
    sense = 0;
  if(sense==2)
    cerr << "C_R_faceDir - error face not adajacent to region\n";
  return sense;
}

/*-------------------------------------------------------------------------
  return sense of face use in region:
  -1 = face is not used by region 
   0 = face normal points into region
   1 = face normal points out of region
   2 = used both ways, face embedded in region
-------------------------------------------------------------------------*/
int C_R_useFace(void *robj, void *fobj) 
{
  GRegion *r = (GRegion*)robj;
  GFace *f = (GFace*)fobj;

  int status = -1;
  if( f->region(1) == r)
    status = 0;
  if( f->region(0) == r)
    status = status == 0 ? 2 : 1;
  return status;
}

/* return pointer to model region on 'dir' side of face 'ftag' */
void *C_F_Region(void *ftag, int dir)
{
  GRegion *region = ((GFace *)ftag)->region(dir);
  return (void*)region;
}

/* return modeler tolerence associated w/ a model object ptr */
int C_modtol(int type, void *gobj, double *tolerance)
{
  *tolerance = ((AOMD::GEntity *)gobj)->tolerance() ;
  return 1 ;
}

int C_getunv(double *length, double *origin)
{
  AOMD::SBoundingBox3d box = theModel->bounds();
  box.makeCube();
  box *= 1.1;
  AOMD::SPoint3 orig(box.min());
  AOMD::SPoint3 len(box.max()-box.min());
  orig.position(origin);
  len.position(length);
  return 1;
}

int C_chkContainment(int type, void *ent, int iflag, double *pnt)
{
  switch(iflag){
  case 1:{ // global coord given
    AOMD::SPoint3 pt(pnt);
    return ((AOMD::GEntity*)ent)->classifyPoint(pt,1,nullptr);
    break;
  }
  case 2: // parametric coord given
    switch(type){
    case 1:
      return ((GEdge*)ent)->classifyParam(*pnt,1,nullptr);
      break;
    case 2:{
      SPoint2 par(pnt);
      return ((GFace*)ent)->classifyParam(par,1,nullptr);
      break;
    }
    default:
      cerr << "C_chkContainment - input error\n";
      break;
    }
  default:
    cerr << "C_chkContainment - input error\n";
    break;
  }
  return 0;
}

int C_chkPtOnBoundary(int type, void *ent, int iflag, double *pnt, 
                     int *bndtyp, void **bndent)
{
  switch(iflag){
  case 1:{ // global coord given
    AOMD::SPoint3 pt(pnt);
    return ((AOMD::GEntity*)ent)->classifyPoint(pt,1,(AOMD::GEntity**)bndent);
    break;
  }
  case 2: // parametric coord given
    switch(type){
    case 1:
      return ((GEdge*)ent)->classifyParam(*pnt,1,(AOMD::GEntity**)bndent);
      break;
    case 2:{
      SPoint2 par(pnt);
      return ((GFace*)ent)->classifyParam(par,1,(AOMD::GEntity**)bndent);
      break;
    }
    default:
      cerr << "C_chkContainment - input error\n";
      break;
    } 
  default:
    cerr << "C_chkContainment - input error\n";
    break;
  }
  return 0;
}

int C_param(int type, void *tag, double *xyz, double *par) 
{
  AOMD::SPoint3 pt(xyz);
  switch(type){
  case 1:{
    *par = ((GEdge*)tag)->parFromPoint(pt);
    if(*par==MAXDOUBLE)
      return 0;
    return 1;
  }
  case 2:{
    SPoint2 ptpar = ((GFace*)tag)->parFromPoint(pt);
    ptpar.position(par);
    if(par[0]==MAXDOUBLE)
      return 0;
    return 1;
  }
  default:
    cerr << "Error - C_param called for region or vertex\n";
  }
  return 0;
}

int C_parType2(int type, void *tag, int dim,int *nbr, double par[2])
{
  if(type!=2){
    cerr << "C_parType2 called for other than face \n";
    return 0;
  }
  GFace *f = (GFace*)tag;
  *nbr = f->paramDegeneracies(dim,par);
  return 1;
}

#ifdef __cplusplus
}  // bracket for extern "C"
#endif
