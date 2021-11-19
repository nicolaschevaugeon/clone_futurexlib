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
#ifndef SIM
#include <iostream>
#include "AOMD.h"
#include "AOMDInternals.h"
#include "mEntity.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mFace.h"
#include "mRegion.h"
#include "mPoint.h"
#include "mException.h"
#include "mAttachableDataContainer.h"
#include "mAOMD.h"
#include "mMesh.h"
// #include "LagrangeMapping.h"
#include "mIteratorNew.h"
#include "ParUtil.h"
#include "PList.h"
#include "AOMDList.h"
#include "mBuildAdj.h"
#include "modeler.h"

#include <cstdio>
#include <cassert>
#include <list>
#include <cstring>
#include "AOMD_cint.h"

using namespace AOMD;
using std::cout;


/*attachable data's*/
class mAttachableVoid : public mAttachableData
{
  public :
    ~mAttachableVoid () override
    {
      //      printf("possible memory leak in a data attachement \n");
      //      delete veryNastyPointer;
    }
    void *veryNastyPointer;
};

pMesh MS_newMesh(pGModel pm)
{
  mMesh *m = new mMesh(1,pm);
  //  allMeshes.push_back(m);  
  return m;
}

pGModel M_model(pMesh mesh)
{
  return mesh->getSolidModel();
}

void P_setParam1(pPoint, double param)
{
  throw 1;
}
void P_setParam2(pPoint, double p1, double p2, int ptch)
{
  throw 1;
}

void M_delete (pMesh pm)
{
//  delete pm;
}

pPList  AOMD_createAdjPList ( pEntity pe , int dim )
{
  pPList pl = PList_new();
  if(pe->getLevel() >= dim || pe->isAdjacencyCreated(dim))
    {  
      for(int i=0;i<pe->size(dim);i++)
	PList_append(pl,(void *)pe->get(dim,i));
    }
  else
     {
      mAdjacencyContainer upward;
      pe->getHigherOrderUpward (dim,upward);
      for(int i=0;i<upward.size();i++)
	PList_append(pl,(void *)upward[i]);
    }
  return pl;
}

//  pPList fakeListToPList(const fakeList *fl)
//  {
//    pPList pl = PList_new();
//    for( fakeList::iter it=fl->begin(); it!=fl->end(); ++it  )
//      PList_append(pl,(void *)(*it));
//    return pl;
//  }
// int M_load_old (pMesh pm, const char *filename)
// {
// #ifdef PARALLEL
//   ParUtil::Instance()->Msg(ParUtil::INFO,"\n***************\n  old PAOMD M_load\n");  
//   myZoltan_LoadBalancer zlb;
//   AOMD_Util::Instance()->import_and_loadbalance(filename,pm,zlb);
// #else
//   cout<<"\n***************\n  old AOMD M_load\n";
//   AOMD_Util::Instance()->import(filename,pm);  
//   M_SetRepresentationToOneLevel(pm);
//   pm->printNumEntities();
// #endif
//   /// into import -> build a null model if necessary
//   return 1;
// }


int M_load (pMesh pm, const char *filename)
{
#ifdef DMUM
// pause DMUM temporarily
  bool wasOn=false;
  if (pm->DMUM_on) 
  {
    wasOn=true;
    pm->DMUM_on=false;
  }
#endif  
  double t1 = ParUtil::Instance()->wTime();
  if (ParUtil::Instance()->rank()==0)
    printf("\n***************\n  AOMD v1_6_0\n");  
    
  char ext[6];
  strcpy(ext,filename+(strlen(filename)-4));
  if(!strcmp(ext,".sms") || !strcmp(ext,".aom")) // M_load
  {

    AOMD_Util::Instance()->import_oneLevel(filename,pm);    

  }
  else  // old M_load done by JF
  {
//     M_load_old(pm, filename);
      throw;
  }  
  double t2 = ParUtil::Instance()->wTime();
  if (ParUtil::Instance()->rank()==0)
  {  
    printf("   (t = %f sec)\n", t2-t1);
    printf("***************\n");  
  }
#ifdef DMUM
  if (wasOn) pm->DMUM_on=true;
#endif  

  return 1;
}

void M_save_topology(char *fname, pMesh m)
{
  AOMD_Util::Instance()->print_topology(fname,m);
}


void M_writeSMS(pMesh pm, const char *name, int version)
{
  char realName [256];
  char ext[1];
  strcpy(ext,name+(strlen(name)-4));
  if(!strcmp(ext,"."))
  {

  sprintf(realName,"%s.sms",name);

  AOMD_Util::Instance()->ex_port(realName,pm,true);
  }
  else AOMD_Util::Instance()->ex_port(name,pm,true);
}

int M_numRegions (pMesh pm)
{
  return pm->size(3);
}

int M_numFaces (pMesh pm)
{
  return pm->size(2);
}

int M_numEdges (pMesh pm)
{
  return pm->size(1);
}

int M_numVertices (pMesh pm)
{
  return pm->size(0);
}

int M_numPoints(pMesh)
{
  throw 1;
}

pRegion M_region(pMesh pm, int i)
{
  if(!i)return *(pm->begin(3));
  else throw new mException(__LINE__,__FILE__,"No random access iterators to regions in AOMD");
}

pFace M_face(pMesh pm, int i)
{
  if(!i)return *(pm->begin(2));
  else throw new mException(__LINE__,__FILE__,"No random access iterators to faces in AOMD");
}

pRegion M_edge(pMesh pm, int i)
{
  if(!i)return *(pm->begin(1));
  else throw new mException(__LINE__,__FILE__,"No random access iterators to edges in AOMD");
}

pVertex M_vertex(pMesh pm, int i)
{
  if(!i)return *(pm->begin(0));
  else throw new mException(__LINE__,__FILE__,"No random access iterators to vertices in AOMD");
}

RIter M_regionIter(pMesh mesh)
{
  return mesh->getEntities().NewIterator (3);
}

FIter M_faceIter(pMesh mesh)
{
  return mesh->getEntities().NewIterator (2);
}

EIter M_edgeIter(pMesh mesh)
{
  return mesh->getEntities().NewIterator (1);
}

VIter M_vertexIter(pMesh mesh)
{
  return mesh->getEntities().NewIterator (0);
}

pRegion RIter_next(RIter it)
{
  if(it->end())return nullptr;
  mEntity *e = *(*it);
  it->next();
  return e;
}

void RIter_delete(RIter it)
{
  it->getEntities()->DeleteIterator(it);
}

void RIter_reset(RIter it)
{
  it->reset();
}

pFace FIter_next(FIter it)
{
  return RIter_next(it);
}

void FIter_delete(FIter it)
{
  RIter_delete(it);
}

void FIter_reset(FIter it)
{
  RIter_reset(it);
}

pEdge EIter_next(EIter it)
{
  return RIter_next(it);
}

void EIter_delete(EIter it)
{
  RIter_delete(it);
}

void EIter_reset(EIter it)
{
  RIter_reset(it);
}

pVertex VIter_next(VIter it)
{
  return RIter_next(it);
}

void VIter_delete(VIter it)
{
  RIter_delete(it);
}

void VIter_reset(VIter it)
{
  FIter_reset(it);
}

/*
 Regions operators
*/

int R_numFaces(pRegion pr)
{
  assert (pr->getLevel() == 3);
  return pr->size(2);
}

pFace R_face(pRegion pr, int n)
{
  return pr->get(2,n);
}

int R_numEdges(pRegion pr)
{
  return pr->size(1);
}

pEdge R_edge(pRegion pr, int n)
{
  return pr->get(1,n);
}

int R_numVertices(pRegion pr)
{
  return pr->size(0);
}

pVertex R_vertex(pRegion pr, int n)
{
  return pr->get(0,n);
}

int R_faceDir (pRegion pr, int n)
{
  return pr->getUse ( pr->get (2,n) );
}

int F_edgeDir (pFace pf, int n)
{
  return pf->getUse (pf->get(1,n));
}

pGRegion R_whatIn(pRegion pr)
{
  if (pr->getClassification()) 
    return (pGRegion)pr->getClassification();
  else
    return (pGRegion)nullptr;
}

pGEntity EN_whatIn(pEntity pe)
{
  return (pGEntity)R_whatIn (pe);
}

pGEntity F_whatIn(pFace pf)
{
  return (pGEntity)R_whatIn (pf);
}

pGEntity  E_whatIn(pEdge pe)
{
  return (pGEntity)R_whatIn (pe);
}

pGEntity  V_whatIn(pVertex pv)
{ 
  return (pGEntity)R_whatIn (pv);
}


int F_numEdges(pFace pf)
{
  return pf->size(1);
}

pEdge F_edge(pFace pg, int n)
{
  return pg->get(1,n);
}

int F_numVertices(pFace pf)
{
  return pf->size(0);
}

pVertex F_vertex(pFace pf, int n)
{
  return pf->get(0,n);
}

int E_numRegions(pEdge pv)
{
  AOMD::mAdjacencyContainer upward;
  pv->getHigherOrderUpward(3,upward);
  return upward.size();
}

int E_numVertices(pEdge pe)
{
  return pe->size(0);
}

pVertex E_vertex(pEdge pe, int n)
{
  //  printf("asking vertex %d of an edge ",n);pe->print();
  return pe->get(0,n);
}

pPList E_vertices(pEdge pe)
{
  pPList pl = PList_new();
  PList_append(pl,(void *)pe->get(0,0));
  PList_append(pl,(void *)pe->get(0,1));
  return pl;
}

int EN_type(pEntity pe)
{
  return pe->getLevel();
}

pFace E_face(pEdge pe, int n)
{
  return pe->get(2,n);
}

int E_numFaces(pEdge pe)
{
  if(pe->isAdjacencyCreated(2))
    return pe->size(2);
  else
    return 0;
}

RIter M_classifiedIter(pMesh mesh, pGEntity ent, int dim)
{
  /// look here for how to get dimension and id of the pGEntity
  /// in order to retrieve gEntity from pEntity
  
//  GEntity *gent = mesh->getGEntity (GEN_tag(ent),GEN_type(ent)); /* GEN_id 
//is changed to GEN_tag since tag is the identification for the gEntities in
//AOMD now - Mohan Nuggehally 06-12-03 */

//  assert (gent->getSolidModelEntity() == ent);
  return mesh->getEntities().NewIterator (dim,ent);
  return nullptr;
}

int F_numRegions(pFace pf)
{
  if(pf->isAdjacencyCreated(3))
    return pf->size(3);
  else
    return 0;
}

RIter M_classifiedRegionIter(pMesh mesh, pGEntity ent)
{
  return M_classifiedIter(mesh,ent, 3);
}

FIter M_classifiedFaceIter(pMesh mesh, pGEntity ent, int closure)
{
  if(!closure)
    return M_classifiedIter(mesh,ent, 2);
  else
    throw 1;
}

EIter M_classifiedEdgeIter(pMesh mesh, pGEntity ent, int closure)
{
  if(!closure)
    return M_classifiedIter(mesh,ent, 1);
  else
    throw 1;
}

VIter M_classifiedVertexIter(pMesh mesh, pGEntity ent, int closure)
{
  if(!closure)
    return M_classifiedIter(mesh,ent, 0);
  else
    throw 1;   
}


pPList F_vertices(pFace pf, int dir)
{
  if (dir <= 0)
    {
      pPList pl = PList_new();
      int N = pf->size(0); 
      for (int i=0;i<N;i++)
	{
	  PList_append(pl,pf->get(0,N-i-1));
	}
      return pl;
    }
  else
    return AOMD_createAdjPList ( pf, 0 );
}

pPList V_regions(pVertex v)
{
  return AOMD_createAdjPList ( v, 3 );
}

pPList E_regions(pEdge e)
{
  return AOMD_createAdjPList ( e, 3 );
}

pPList V_faces(pVertex v)
{
  return AOMD_createAdjPList ( v, 2 );
}


int EN_id(pEntity pe)
{
  if(pe->getLevel() == 0)return pe->getId();
  else return pe->getAttachedInt(AOMD_Util::Instance()->getId());
}

void EN_setID(pEntity pe, int id)
{
  if (pe->getLevel() == 0)
  throw new mException(__LINE__,__FILE__,"Can't set ID to vertices in AOMD");
  pe->attachInt(AOMD_Util::Instance()->getId(),id);
}

pPList R_edges(pRegion pr,int)
{
  return AOMD_createAdjPList ( pr, 1 );
}


pPList R_faces(pFace pr, int)
{
  return AOMD_createAdjPList ( pr, 2 );
}

pPList R_vertices(pRegion pr, int)
{
  return AOMD_createAdjPList ( pr, 0 );
}

pPoint  V_point(pVertex pv)
{
  return pv;
}

pPList F_edges(pFace pf , int dir,pVertex vtx)
{
  return AOMD_createAdjPList ( pf, 1 );
}

pPList F_regions(pFace pf)
{
  return AOMD_createAdjPList ( pf, 3 );
}

pRegion F_region(pFace pf, int n)
{
  if(!pf->isAdjacencyCreated(3))return nullptr;
  if(n >= pf->size(3))return nullptr;
  if(pf->size(3) > 2)throw;
  return pf->get(3,n);
}

int V_numEdges(pVertex  pv)
{
  if(pv->isAdjacencyCreated(1))
    return pv->size(1);
  else
    return 0;
}

pEdge V_edge(pVertex pv, int n)
{
  return pv->get(1,n);
}

pPList V_edges(pVertex pv)
{
  return AOMD_createAdjPList ( pv, 1 );
}

int V_numFaces(pVertex pv)
{
  AOMD::mAdjacencyContainer upward;
  pv->getHigherOrderUpward(2,upward);
  return upward.size();
}

int V_numRegions(pVertex pv)
{
  AOMD::mAdjacencyContainer upward;
  pv->getHigherOrderUpward(3,upward);
  return upward.size();
}

double P_x(pPoint p){return ((mVertex*)p)->point()(0);}
double P_y(pPoint p){return ((mVertex*)p)->point()(1);}
double P_z(pPoint p){return ((mVertex*)p)->point()(2);}

pMeshDataId MD_newMeshDataId(const char *tag) {return AOMD_Util::Instance()->newMeshDataId(tag);}
pMeshDataId MD_lookupMeshDataId(const char *tag) {return AOMD_Util::Instance()->lookupMeshDataId   (tag);}
void MD_deleteMeshDataId(pMeshDataId id) {AOMD_Util::Instance()->deleteMeshDataId(id);}

void EN_attachDataI(pEntity pe, const char *tag, int data)
{
  pe->attachInt(MD_lookupMeshDataId(tag),data);
}

void EN_attachDataInt(pEntity pe, pMeshDataId id, int data)
{
  pe->attachInt(id,data);
}

int EN_dataI(pEntity pe, const char *tag)
{
  return pe->getAttachedInt(MD_lookupMeshDataId(tag));
}

int EN_getDataDbl(pEntity ent, pMeshDataId id, double *value)
{
  *value = ent->getAttachedDouble(id);
  return 1;
}

void EN_attachDataDbl(pEntity ent, pMeshDataId id, double value)
{
  ent->attachDouble(id,value);
}

int EN_modifyDataI(pEntity pe, const char *tag, int data)
{
  pe->attachInt(MD_lookupMeshDataId(tag),data);
  return 1;
}

void EN_modifyDataInt(pEntity ent, pMeshDataId id, int value)
{
  ent->attachInt(id,value);
}

void MS_init(void)
{
}

void MS_exit(void)
{
}

double P_param1(pPoint p)
{
  if (!p->getData(AOMD_Util::Instance()->getParametric()))
    throw new mException(__LINE__,__FILE__,
			 "The mesh you have loaded does not have parametric coorinates");
  Trellis_Util::mVector vec = p->getAttachedVector (AOMD_Util::Instance()->getParametric());
  return vec(0);
}
void P_param2(pPoint p, double *p1, double *p2, int *ptch)
{
  if (!p->getData(AOMD_Util::Instance()->getParametric()))
    throw new mException(__LINE__,__FILE__,
			 "The mesh you have loaded does not have parametric coorinates");
  Trellis_Util::mVector vec = p->getAttachedVector (AOMD_Util::Instance()->getParametric());
  *p1   = vec(0);
  *p2   = vec(1);
//  *ptch = (int)vec(2);
}

pVertex  E_otherVertex(pEdge pe, pVertex v)
{
  mEdge *e = (mEdge*)pe;
  if(e->vertex(0) == v)return e->vertex(1);
  if(e->vertex(1) == v)return e->vertex(0);
  throw new mException(__LINE__,__FILE__,
		       "E_otherVertex failed");
}

pFace E_otherFace(pEdge edge, pFace face, pRegion region)
{
  mEdge *e = (mEdge*)edge;
  mFace *f = (mFace*)face;
  mRegion *r = (mRegion*)region;

  for(int i=0;i<r->size(2);i++)
    {
      mFace *f2 = (mFace*)r->get(2,i);
      if(f2 != f && f2->find(edge))return f2;
    }  
  throw new mException(__LINE__,__FILE__,
		       "E_otherFace failed");
}

void V_coord(pVertex pv, double *xyz)
{
   Trellis_Util::mPoint p = ((mVertex*)pv)->point();
   xyz[0] = p(0);
   xyz[1] = p(1);
   xyz[2] = p(2);
}


pPoint P_new(void)
{
  return new mVertex(0,Trellis_Util::mPoint(0,0,0),nullptr);
}

void P_delete(pPoint p)
{
  delete p;
}

void P_setPos(pPoint p , double x, double y, double z)
{
  Trellis_Util::mPoint temp(x,y,z);
  ((mVertex*)p)->move(Trellis_Util::mPoint(x,y,z));
}

void AOMD_P_setParametricPos ( pPoint p, double p1, double p2, double p3)
{
  p->attachVector (AOMD_Util::Instance()->getParametric() , Trellis_Util::mVector(p1,p2,p3));
}

void EN_attachDataP(pEntity pe, const char *tag, void *data)
{
  unsigned int itag = MD_lookupMeshDataId(tag);
  EN_attachDataPtr(pe, itag, data);
}

void EN_attachDataPtr(pEntity ent, pMeshDataId id, void * value)
{
  mAttachableVoid *av = (mAttachableVoid *)ent->getData(id);
  if(!av)
    {
      av = new mAttachableVoid;
      ent->attachData(id,av);
    }
  else
    {
      //      delete av->veryNastyPointer;
    }
  //  printf ("attaching %s to ",tag);pe->print();
  av->veryNastyPointer = value;
}

void * EN_dataP(pEntity pe, const char *tag)
{  
  unsigned int itag = MD_lookupMeshDataId(tag);
  mAttachableVoid *av = (mAttachableVoid *)pe->getData(itag);
  if(!av)
    {
      //      printf("data %s not found on ",tag);pe->print();
      //      pe->printAllAttachable ();
    }
  if(!av)return nullptr;
  return av->veryNastyPointer;
}

int EN_modifyDataP(pEntity pe, const char *tag, void * data)
{
  EN_attachDataP(pe, tag, data);
  return 1;
}

// void F_normalVector(pFace face, int dir, double* normal)
// {
//   Trellis_Util::LagrangeMapping mapping(face);
//   double u,v,w;
//   Trellis_Util::mTensor2 j;
//   mapping.COG(u,v,w);
//   mapping.jacInverse(u,v,w, j);
//   
//   // in reference coordinates, n = (0,0,1);
//   Trellis_Util::mVector n (0,0,1);
//   n*= j;
//   
//   normal[0] = n(0) * dir;
//   normal[1] = n(1) * dir;
//   normal[2] = n(2) * dir;  
// }

int E_numPoints(pEdge)
{
  return 0;
}

pPoint E_point(pEdge , int n)
{
  throw new mException (__LINE__, __FILE__, "points are not supported on mesh edges, ask jf");  
}

int EN_getDataPtr(pEntity ent, pMeshDataId id, void **value)
{

  mAttachableVoid *av = (mAttachableVoid *)ent->getData(id);
  if(!av)
    {
      //      printf("data %s not found on ",tag);pe->print();
      //      pe->printAllAttachable ();
    }
  if(!av)return 0;
  *value =  av->veryNastyPointer;
  return 1;
}

int EN_getDataInt(pEntity ent, pMeshDataId id, int *value)
{
  *value = ent->getAttachedInt(id);
  if(*value)return 1;
  return 0;
}

void R_setWhatIn(pRegion e, pGEntity what)
{
      e->classify (what);
      return;
}

void EN_setWhatIn(pMesh pm, pEntity e, pGEntity what)
{
  e->classify (what);
}

void F_setWhatIn(pFace   e, pGEntity what){R_setWhatIn(e,what);}
void E_setWhatIn(pEdge   e, pGEntity what){R_setWhatIn(e,what);}
void V_setWhatIn(pVertex e, pGEntity what){R_setWhatIn(e,what);}

pRegion M_createR(pMesh mesh, int nFace, pFace *faces, int *dirs, pGEntity gent)
{
  if(nFace == 4) 
    {
      pRegion r = (mRegion*)mesh->createTetWithFaces ((mFace*)faces[0],(mFace*)faces[1],
				(mFace*)faces[2],(mFace*)faces[3],gent);
      faces[0]->add(r);
      faces[1]->add(r);
      faces[2]->add(r);
      faces[3]->add(r);
#ifdef DEBUG
      if (R_Volume((mEntity*)r)<0.0)
        std::cerr<<"AOMD FATAL: negative volume region created "<<r->getUid()<<std::endl;
      assert (R_Volume((mEntity*)r)>0.0);
#endif      
      return r;
    }
  else return nullptr;
}

pFace M_createF(pMesh mesh, int nEdge, pEdge *edges, int *dirs, pGEntity gent)
{
  pFace f;
  if(nEdge == 3) 
    f = mesh->createFaceWithEdges((mEdge*)edges[0],(mEdge*)edges[1],
 			          (mEdge*)edges[2],gent,dirs);
  else if (nEdge == 4)
    f =  mesh->createFaceWithEdges((mEdge*)edges[0],(mEdge*)edges[1],
			           (mEdge*)edges[2],(mEdge*)edges[3],gent,dirs);
  else 
    return nullptr;
  for(int i=0;i<nEdge;i++)
    edges[i]->add(f);
  return f;
}

pEdge M_createE(pMesh mesh, pVertex v1, pVertex v2, pGEntity gent)
{
#ifdef DEBUG
  pEdge e = E_exist(v1, v2);
  if (e)
  {
    std::cerr<<"AOMD::FATAL: Attemp to create existing edge "<<e->getUid()<<"\n  => ";
    e->print();
    assert(!e);
  }       
#endif
  pEdge p = (pEdge)mesh->createEdge ((mVertex*)v1,(mVertex*)v2,gent);
  v1->add (p);
  v2->add (p);
  return p;
}

pVertex M_createVP(pMesh mesh, double x, double y, double z, double *param,
		   int id, pGEntity gent)
{
  pVertex vv;
  if (id != 0)
    vv=(pVertex)mesh->createVertex (id,x,y,z,gent);
  else
    vv=(pVertex)mesh->createVertex (x,y,z,gent);

  vv->attachVector ( AOMD_Util::Instance()->getParametric() , Trellis_Util::mVector (*param,*(param+1),*(param+2)) );
  return vv;
}

pVertex M_createVP2(pMesh mesh, double *xyz, double *param,
		    int id, pGEntity gent)
{
  pVertex vv;
  if (id!=0)
    vv=mesh->createVertex (id,xyz[0],xyz[1],xyz[2],gent);
  else 
    vv=mesh->createVertex (xyz[0],xyz[1],xyz[2],gent);
  vv->attachVector ( AOMD_Util::Instance()->getParametric() , Trellis_Util::mVector (*param,*(param+1),*(param+2)) );
  return vv;
}

pPoint M_createP(pMesh mesh, double x, double y, double z, double *param,
		 int unused, pGEntity ent)
{
  throw 1;
}

void M_remove(pMesh m, pRegion region, bool do_delete = true) 
{
  int dim = region->getLevel();  
  if(region->isAdjacencyCreated(dim-1))
    {
      for(int i=0;i<region->size(dim-1);i++)
	region->get(dim-1,i)->del(region);  
    }
  if (do_delete)m->DEL(region);
}

void M_removeRegion(pMesh m, pRegion r){M_remove(m,r);}

void M_removeFace(pMesh m, pFace face)
{
//#ifdef PARALLEL
//  if (m->size(3)>0 && face->getCommonBdry())
//    EN_deleteRemoteCopy(m,face);
//#endif
M_remove(m,face);
}

void M_removeEdge(pMesh m, pEdge edge)
{
//#ifdef PARALLEL
//  if (edge->getCommonBdry())
//    EN_deleteRemoteCopy(m,edge);
//#endif
  M_remove(m,edge);
}

void M_removeVertex(pMesh m, pVertex vertex)
{
//#ifdef PARALLEL
//  if (vertex->getCommonBdry())
//    EN_deleteRemoteCopy(m,vertex);
//#endif
//  cout<<"*** ("<<P_pid()<<") Remove "<<vertex->getUid()<<"\n";
  m->DEL(vertex);
}

void M_removePoint(pMesh, pPoint point){}


void F_chDir(pFace pf) 
{
  pVertex v0 = pf->get(0,0);
  pVertex v1 = pf->get(0,1);
  pVertex v2 = pf->get(0,2);
  mEdge* e0 = (mEdge*)pf->get(1,0);
  mEdge* e1 = (mEdge*)pf->get(1,1);
  mEdge* e2 = (mEdge*)pf->get(1,2);
  pf->deleteAdjacencies(0);
  pf->deleteAdjacencies(1);
  pf->add(e0);
  pf->add(e2);
  pf->add(e1);
  pf->add(e0->commonVertex(e1));
  pf->add(e0->commonVertex(e2));
  pf->add(e1->commonVertex(e2));
}

int EN_whatInType(pEntity e)
{
  /// do something here
  return GEN_type(e->getClassification()); 
}
int R_whatInType(pRegion e) { return EN_whatInType(e);}
int F_whatInType(pFace e){ return EN_whatInType(e);}
int E_whatInType(pEdge e){ return EN_whatInType(e);}
int V_whatInType(pVertex e){ return EN_whatInType(e);}

void EN_removeData(pEntity pe, const char *tag)
{
  unsigned int itag = MD_lookupMeshDataId(tag);
  pe->deleteData (itag);
}

void EN_deleteData(pEntity pe, pMeshDataId id)
{
  pe->deleteData (id);
}

int E_inClosure(pEdge pe, pEntity ent)
{
  return R_inClosure(pe, ent);
}

int R_dirUsingFace(pRegion pr, pFace face)
{
  return pr->getUse(face);
}

int F_dirUsingEdge(pFace pf, pEdge edge)
{
   return pf->getUse(edge);
}

void EN_modifyDataPtr(pEntity ent, pMeshDataId id, void * value)
{
  EN_attachDataPtr(ent, id, value);
}


pRegion M_nextRegion(pMesh, void **restart)
{
  throw 1;
}

int M_numClassifiedEntities(pMesh mesh, pGEntity ent, int dim)
{
  mMesh::iter it = mesh->begin(dim);
  mMesh::iter ite = mesh->end(dim);
  int k=0;
  while(it != ite)
    {
      //if((*it)->getClassification()) k++; 
      if((*it)->getClassification() == ent) k++;
      ++it;
    }
  return k;
}

int M_numClassifiedRegions(pMesh mesh, pGEntity ent)
{
  return M_numClassifiedEntities(mesh, ent, 3);
}
int M_numClassifiedFaces(pMesh mesh, pGEntity ent)
{
  return M_numClassifiedEntities(mesh, ent, 2);
}
int M_numClassifiedEdges(pMesh mesh, pGEntity ent)
{
  return M_numClassifiedEntities(mesh, ent, 1);
}
int M_numClassifiedVertices(pMesh mesh, pGEntity ent)
{
  return M_numClassifiedEntities(mesh, ent, 0);
}

void V_setSize(pVertex pe,double s)
{
  pe->attachDouble(AOMD_Util::Instance()->getSize(),s);
}

double V_size (pVertex pe)
{
  return pe->getAttachedDouble(AOMD_Util::Instance()->getSize());
}

static unsigned int markIdsFirst = 1000000000;

MarkID MD_registerMark(char *name, unsigned int size, unsigned int flags)
{
  return 1;
}

MarkID MD_markID(char *name)
{
  return 1;
}
unsigned int EN_mark(pEntity entity, MarkID id)
{
  return entity->getAttachedInt(markIdsFirst+id);  
}

void EN_setMark(pEntity entity, MarkID id, unsigned int value)
{
  entity->attachInt(markIdsFirst+id,value);  
}

void V_setPoint(pVertex v, pPoint pt)
{
  ((AOMD::mVertex*)v)->move (((AOMD::mVertex*)pt)->point());
}

int R_inClosure(pRegion pr, pEntity ent)
{
  fakeList f (pr, ent->getLevel());
  return  f.inlist (ent);
}

int F_inClosure(pFace f, pEntity ent)
{
  return R_inClosure(f, ent);
}

void V_merge2(pMesh pM, pVertex, pVertex)
{
  throw 1;
}
void V_merge(pVertex, pVertex)
{
  throw 1;
}

void E_merge2(pMesh pM, pEdge, pEdge)
{
  throw 1;
}
void E_merge(pEdge, pEdge)
{
  throw 1;
}
void F_merge2(pMesh pM, pFace, pFace)
{
  throw 1;
}
void F_merge(pFace, pFace)
{
  throw 1;
}

struct four_ints
{
  int v[4];
  four_ints () = default;;  
  four_ints (pEntity e)    
  {
    v[0] = v[1] = v[2] = v[3] = -1;
    for (int i=0;i<e->size(0);i++)
      {
	v[i] = e->get(0,i)->getId();
      }
    std::sort(v,v+4);
  }
  four_ints (pEntity e, pEntity pV, pEntity pV2)    
  {
    v[0] = v[1] = v[2] = v[3] = -1;
    for (int i=0;i<e->size(0);i++)
      {
	v[i] = (e->get(0,i) == pV)?pV2->getId():e->get(0,i)->getId();
      }
    std::sort(v,v+4);
  }
};

struct lessthan_four_ints
{
  bool operator () (const four_ints &v1, const four_ints &v2) const
  {
    if (v1.v[0] < v2.v[0])return true;
    if (v1.v[0] > v2.v[0])return false;
    if (v1.v[1] < v2.v[1])return true;
    if (v1.v[1] > v2.v[1])return false;
    if (v1.v[2] < v2.v[2])return true;
    if (v1.v[2] > v2.v[2])return false;
    if (v1.v[3] < v2.v[3])return true;
    return false;
  }
};

static pGEntity getClassif (std::map<four_ints,AOMD::GEntity*, lessthan_four_ints> &cr,
			   pEntity e)
{
  std::map<four_ints,AOMD::GEntity*, lessthan_four_ints>::iterator it =
    cr.find(four_ints(e));
  if (it == cr.end())
    return e->getClassification();
  return it->second;
}

static void updateClassif (std::map<four_ints,AOMD::GEntity*, lessthan_four_ints> &cr,
			       pEntity e, pEntity pV1 = nullptr , pEntity pV2  = nullptr)
{
  four_ints xxx;

  if (pV1)
    xxx = four_ints (e,pV1,pV2);
  else
    xxx = four_ints (e);
    
  std::map<four_ints,AOMD::GEntity*, lessthan_four_ints>::iterator it =
    cr.find(xxx);
  if (it == cr.end())
    {
      cr[xxx] = e->getClassification();
    }
  else
    {
      pGEntity ge = (*it).second;
      if (GEN_type(ge) > GEN_type(e->getClassification()))
	cr[xxx] = e->getClassification();
    }
}

pPList E_faces(pEdge e) 
{
  return AOMD_createAdjPList (e, 2);
}

int M_findE(pMesh mesh, pEntity ent, int i)
{
  if (std::find(mesh->beginall(i), mesh->endall(i), ent)!= mesh->endall(i))
    return 1;
  return 0;    
}
#endif
