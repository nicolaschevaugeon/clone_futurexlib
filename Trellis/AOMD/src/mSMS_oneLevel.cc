/**************************************************************************** 
 * 
 *   Copyright (C) 2004
 *   Rensselaer Polytechnic Institute
 * 
 *   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written 
 *   and maintained by the Scientific Computation Research Center (SCOREC) at 
 *   Rensselaer Polytechnic Intitute, Troy, NY, USA.
 * 
 *   This program is free software; you can redistribute it and/or modify it
 *   under the terms of the Rensselaer SCOREC Public License.
 * 
 *   This program is distributed in the hope that it will be useful, 
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   license text for more details.
 *   
 *   You should have received a copy of the Rensselaer SCOREC Public License
 *   along with this program; if not, write to Rensselaer Polytechnic Institure,
 *   110 8th Street, SCOREC, Troy, NY  12180, USA
 * 
 *****************************************************************************/

/*****************************************************************************
 * 
 *  src/mSMS_oneLevel.cc
 *  Created by Eunyoung Seol, on Mon Oct 29 2003, 10:26:37 EDT
 *
 *  File Content: symmetrix mesh import/export finctions with onelevelRepresentation
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdio>
#include "mAOMD.h"
#include "AOMD_Internals.h"
#include "AOMD_cint.h"

#include "mMesh.h"
#include "mVertex.h"
#include "mEdge.h"
#include "mFace.h"
#include "mTet.h"
#include "mException.h"
#include <cassert>
#include <fstream>

#include "AOMD.h"
#include "modeler.h"
#include "SModelMember.h"
#include "SGModel.h"
#ifdef ACIS
#include "AcisFace.h"
#endif

#ifdef FLEXDB
#include "mFlexDB.h"
#endif

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::vector;

namespace AOMD {
    
    const double accuracy=0.000005;
    
    bool isEqual(double x, double y)
    {  
        if (fabs(x)<accuracy && fabs(y)<accuracy) return true;
        if (fabs((x-y)/(0.5*(x+y))) < accuracy) return true;
        return false;
    }
    
    void AOMD_Util::importSmsFile_oneLevel(const char *fName, mMesh *theMesh)
    {
        FILE *in = fopen (fName,"r");
        if(!in)throw new mException (__LINE__,__FILE__,"unable to open file in loadSmsFile");
        
        char line[256];
        int i,j;
        
        int NbRegions=0,NbFaces,NbEdges,NbVertices,NbPoints;
        int GEntityType,GEntityId,EntityNbConnections;
        int Dummy;
        int VertexId1,VertexId2;
        int Edge1,Edge2,Edge3;
        int Face1, Face2, Face3, Face4;
        int nbPts,NbEdgesOnFace,NbFacesOnRegion;
        double x,y,z;
        
        fscanf(in,"%s %d",line,&Dummy);
        typedef GEntity* (*entityBy_FP)(SGModel*,int,int);
        entityBy_FP fp = GM_entityByTag;
        
        if (Dummy==2)
        {  
            fp = GM_entityByID;
            theMesh->setGEntity_FP(GM_entityByTag);
        }
        #ifdef ACIS
        //  cout<<"accuracy = "<<accuracy<<endl;
        if (Dummy==4)
        {
            // adjust the tag of gEntities for acis model
            int numGV, numGE, numGF, numGR, numDAdj, numUAdj, numV;
            int gType, gTag, adjGtag;
            double xyz[3];
            bool found, failed;
            vector<int> adjDE, adjUE;
            pPList evertices, efaces, fedges, fvertices, rfaces;
            pGVertex gvertex;
            pGEdge gedge;
            pGFace gface;
            pGRegion gregion;
            double u_low, v_low, u_high, v_high;
            
            
            //inFile >>numGV >> numGE >> numGF >>numGR;
            fscanf(in,"%d %d %d %d",&numGV,&numGE,
                   &numGF,&numGR);
            // adjust gv_tag
            for(i=0;i<numGV;i++)
            {
                //inFile >> gType>>gTag>>x>>y>>z;     
                fscanf(in,"%d %d %lf %lf %lf",&gType, &gTag, &x, &y, &z);
                assert(gType==0);
                found=false;
                // iterate over Gvertices and find the corresponding one
                GVIter gvit = GM_vertexIter(theMesh->getSolidModel());
                
                while ((gvertex = GVIter_next(gvit)))
                {
                    GV_point(gvertex, xyz);
                    bool flag=true;
                    if (!isEqual(xyz[0],x)) flag=false;
                    if (!isEqual(xyz[1],y)) flag=false;
                    if (!isEqual(xyz[2],z)) flag=false;
                    if (!flag) 
                    {
                        /*	  cout<<"xyz Read ("<<x<<","<<y<<","<<z<<")\n"
                         *	      <<"xyz of GV "<<GEN_tag(gvertex)<<" ("
                         *	      <<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")\n";
                         */	      
                        continue; // go to next gv
                    }
                    if (!found)
                    {
                        #ifdef DEBUG
                        //          cout<<"  -- set GV"<<GEN_tag((pGEntity)gvertex)<<" to ";
                        #endif
                        ((SModelMember*)gvertex)->setTag(gTag);
                        #ifdef DEBUG
                        /*            cout<<"GV"<<GEN_tag((pGEntity)gvertex)<<"\n";
                         *	    cout<<"\txyz Read ("<<x<<","<<y<<","<<z<<")\n"
                         *	      <<"\txyz of GV "<<GEN_tag(gvertex)<<" ("
                         *	      <<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")\n";
                         */
                        #endif
                        ((SModelMember*)gvertex)->setID(gTag-1);
                        found=true;
                    }
                    else
                    {
                        cout<<"AOMD Error: There are multiple model vertices with given coords\n";
                        cout<<"\txyz Read ("<<x<<","<<y<<","<<z<<")\n"
                        <<"\txyz of GV "<<GEN_tag(gvertex)<<" ("
                        <<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")\n";
                        
                        fclose(in);
                        throw 1;
                    }
                } // while
                
                if (!found)
                {
                    GVIter_reset(gvit);
                    gvertex = GVIter_next(gvit);
                    cout<<"AOMD Error: No model vertex with coords ["
                    <<x*(1/GEN_tolerance((pGEntity)gvertex))
                    <<","<<y*(1/GEN_tolerance((pGEntity)gvertex))
                    <<","<<z*(1/GEN_tolerance((pGEntity)gvertex))<<"] found\n";
                    GVIter_delete(gvit);
                    
                    fclose(in);
                    throw 1;
                }
                GVIter_delete(gvit);      
            } // for numGV
            
            // adjust ge_tag
            for(i=0;i<numGE;i++)
            {
                //      inFile>>gType>>gTag>>numDAdj;     
                fscanf(in,"%d %d %d",&gType,&gTag, &numDAdj);
                
                assert(gType==1);
                found=false;
                adjDE.clear();
                adjUE.clear();
                
                for (j=0; j<numDAdj;++j)
                {
                    //inFile>>adjGtag;
                    fscanf(in,"%d",&adjGtag);
                    
                    adjDE.push_back(adjGtag);      
                }
                
                // read the # faces
                //inFile>>numUAdj;
                fscanf(in,"%d",&numUAdj);
                for (j=0; j<numUAdj;++j)
                {
                    // reaf #vertices of face
                    //inFile>>numV;
                    fscanf(in,"%d",&numV);
                    
                    adjUE.push_back(numV);
                    for (int k=0; k<numV;++k)
                    {
                        // read adj vertex
                        //inFile>>adjGtag;
                        fscanf(in,"%d",&adjGtag);
                        adjUE.push_back(adjGtag);      
                    }
                }
                
                // iterate over Gvertices and find the corresponding one
                GEIter geit = GM_edgeIter(theMesh->getSolidModel());
                while ((gedge = GEIter_next(geit)))
                {
                    failed=false;
                    evertices = GE_vertices(gedge);
                    if (numDAdj != PList_size(evertices))
                        continue;
                    for (j=0; j<numDAdj;++j)  
                    {
                        if (GEN_tag((pGEntity)PList_item(evertices,j))!=adjDE[j])
                        {
                            failed=true;
                            break;
                        }
                    } // for
                    
                    if (failed) continue;
                    
                    efaces = GE_faces(gedge);
                    if (numUAdj != PList_size(efaces))
                        continue;
                    
                    int counter=0;
                    
                    for (j=0; j<numUAdj;++j)  
                    {
                        gface = (pGFace)PList_item(efaces, j);
                        fvertices = GF_vertices(gface);
                        numV=adjUE[counter++];
                        for (int k=0; k<numV;++k)
                        {
                            gvertex = (pGVertex)GM_entityByTag(theMesh->getSolidModel(), 0, adjUE[counter++]);
                            if (!PList_inList(fvertices,(void*)gvertex))
                            {
                                failed=true;
                                break;
                            }  
                        }
                    } // for
                    if (failed) continue;	
                    
                    if (!found)
                    {
                        //            cout<<"  -- set GE"<<GEN_tag((pGEntity)gedge)<<" to ";
                        ((SModelMember*)gedge)->setTag(gTag);
                        ((SModelMember*)gedge)->setID(gTag-1);
                        //            cout<<"GE"<<GEN_tag((pGEntity)gedge)<<"\n";
                        found=true;
                    }
                    else
                    {
                        cout<<"AOMD Error: There are multiple model edges with adj.v+f.\n  adjV: ";
                        fclose(in);
                        throw 1;
                    } // !failed
                }
                GEIter_delete(geit);      
                if (!found)
                {
                    cout<<"AOMD Error: No model edge with given adj.v+f. found\n adjV: ";
                    for (int k=0; k<numDAdj;++k)
                        cout<<"GV"<<adjDE[k]<<", ";
                    cout<<endl;
                    
                    fclose(in);
                    throw 1;
                }
            } // for numGE
            
            
            // adjust gf_tag
            for(i=0;i<numGF;i++)
            {
                //inFile>>gType>>gTag>>numDAdj;     
                fscanf(in,"%d %d %d",&gType, &gTag, &numDAdj);
                assert(gType==2);
                found=false;
                adjDE.clear();adjUE.clear();
                for (j=0; j<numDAdj;++j)
                {
                    //inFile>>adjGtag;
                    fscanf(in,"%d",&adjGtag);
                    adjDE.push_back(adjGtag);
                }
                //inFile>>numV;
                fscanf(in,"%d",&numV);
                for (j=0; j<numV;++j)
                {
                    //inFile>>adjGtag;
                    fscanf(in,"%d",&adjGtag);
                    adjUE.push_back(adjGtag);
                }
                /// read parRange values
                //inFile>>u_low>>u_high>>v_low>>v_high;
                fscanf(in,"%lf %lf %lf %lf",&u_low, &u_high, &v_low, &v_high);
                // iterate over Gvertices and find the corresponding one
                GFIter gfit = GM_faceIter(theMesh->getSolidModel());
                while ((gface = GFIter_next(gfit)))
                {
                    failed=false;
                    fedges = GF_edges(gface);
                    if (numDAdj != PList_size(fedges))
                        continue;
                    for (j=0; j<numDAdj;++j)  
                    {
                        if (GEN_tag((pGEntity)PList_item(fedges,j))!=adjDE[j])
                        {
                            failed=true;
                            break;
                        }
                    } // for
                    if (failed) continue;
                    
                    fvertices = GF_vertices(gface);
                    if (numV != PList_size(fvertices))
                        continue;
                    for (j=0; j<numV;++j)  
                    {
                        if (GEN_tag((pGEntity)PList_item(fvertices,j))!=adjUE[j])
                        {
                            failed=true;
                            break;
                        }
                    } // for
                    if (failed) continue;
                    
                    if (!found)
                    {
                        //          cout<<"  -- set GF"<<GEN_tag((pGEntity)gface)<<" to ";
                        ((SModelMember*)gface)->setTag(gTag);
                        ((SModelMember*)gface)->setID(gTag-1);
                        //eric (bug)	  ((AcisFace*)gface)->setParRange(u_low, u_high, v_low, v_high);
                        //          cout<<"GE"<<GEN_tag((pGEntity)gface)<<"\n";
                        found=true;
                    }
                    else
                    {
                        cout<<"AOMD Error: There are multiple model faces with adj.e.\n";
                        fclose(in);
                        throw 1;
                    }
                }
                GFIter_delete(gfit);      
                if (!found)
                {
                    cout<<"AOMD Error: No model face with given adj.e. found\n";
                    fclose(in);
                    throw 1;
                }
            } // for numGF
            
            // adjust gr_tag
            for(i=0;i<numGR;i++)
            {
                //inFile>>gType>>gTag>>numDAdj;     
                fscanf(in,"%d %d %d",&gType, &gTag, &numDAdj);
                assert(gType==3);
                found=false;
                adjDE.clear();
                for (j=0; j<numDAdj;++j)
                {
                    //inFile>>adjGtag;
                    fscanf(in,"%d",&adjGtag);
                    adjDE.push_back(adjGtag);      
                }
                
                // iterate over Gvertices and find the corresponding one
                GRIter grit = GM_regionIter(theMesh->getSolidModel());
                while ((gregion = GRIter_next(grit)))
                {
                    failed=false;
                    rfaces = GR_faces(gregion);
                    if (numDAdj != PList_size(rfaces))
                        continue;
                    for (j=0; j<numDAdj;++j)  
                    {
                        // let's ignore the order of downward adj.
                        //cout<<"  -- adjV"<<adjE[i]<<endl;
                        if (!PList_inList(rfaces,
                            (void*)GM_entityByTag(theMesh->getSolidModel(), 2, adjDE[j])))
                        {
                            //cout<<"GE"<<GEN_tag((pGEntity)gedge)<<" doesn't match"<<endl;
                            failed=true;
                            break;
                        }
                    } // for
                    if (failed) continue;
                    else
                    {
                        if (!found)
                        {
                            //            cout<<"  -- set GR"<<GEN_tag((pGEntity)gregion)<<" to ";
                            ((SModelMember*)gregion)->setTag(gTag);
                            ((SModelMember*)gregion)->setID(gTag-1);
                            //            cout<<"GR"<<GEN_tag((pGEntity)gregion)<<"\n";
                            found=true;
                        }
                        else
                        {
                            cout<<"AOMD Error: There are multiple model regions with adj.f.\n";
                            fclose(in);
                            throw 1;
                        }
                    } // !failed
                }
                GRIter_delete(grit);      
                if (!found)
                {
                    cout<<"AOMD Error: No model region with given adj.f. found\n";
                    fclose(in);
                    throw 1;
                }
            } // for numGE
            //inFile>>line>>Dummy;
            fscanf(in,"%s %d",line, &Dummy);
            
            if (Dummy==2)
            {  
                fp = GM_entityByID;
                theMesh->setGEntity_FP(GM_entityByTag);
            }
            PList_delete(evertices);
            PList_delete(efaces);
            PList_delete(fedges);
            PList_delete(fvertices);
            PList_delete(rfaces);
        }// if (Dummy==4)
        #endif
        
        
        //  inFile>>NbRegions>>NbFaces>>NbEdges>>NbVertices>>NbPoints;
        fscanf(in,"%d %d %d %d %d",&NbRegions,&NbFaces,&NbEdges,
               &NbVertices,&NbPoints);
        
        #ifdef FLEXDB
        if (NbRegions>0) 
            set_MRM_oneLevel(theMesh,3);
        else 
            set_MRM_oneLevel(theMesh,2);
        theMesh->setMeshModFunctors(true); // one-level
        #endif
        
        if (P_pid()!=0){
            fclose(in);
            return;
        }
        
        cout<<"reading "<<NbVertices<<" vertices, ";
        double u,v;
        int patch;
        
        // *******************************
        //	Read Vertices
        // *******************************
        std::vector<mVertex*> vertices;
        mVertex* vv;        
        for(i=0;i<NbVertices;i++)
        {
            //inFile>>GEntityId;       
            fscanf(in,"%d",&GEntityId); 
            if (GEntityId)
            {
                //inFile>>GEntityType>>EntityNbConnections>>x>>y>>z; 
                fscanf(in,"%d %d %lf %lf %lf",&GEntityType,&EntityNbConnections,
                       &x,&y,&z);
                vv = theMesh->createVertex(i+1,x,y,z,
                                           theMesh->getGEntity(GEntityId-1,GEntityType,fp)); 
                vertices.push_back(vv);
                
                switch(GEntityType)
                {
                    case 1:
                        //inFile>>u;
                        fscanf(in,"%le",&u);
                        vv->attachVector ( getParametric() , Trellis_Util::mVector (u,0,0) );
                        break;
                    case 2:
                        //inFile>>u>>v>>patch;		
                        fscanf(in,"%le %le %d",&u,&v,&patch);	
                        vv->attachVector ( getParametric() , Trellis_Util::mVector (u,v,patch) );
                        break;
                    default: break;
                }
            }
        } // end of reading vertices
        
        #if defined(DMUM) || defined(FLEXDB)
        vv->g_mesh = theMesh;
        #endif 
        // *******************************
        //	Read Edges
        // *******************************
        
        std::vector<mEntity*> edges;
        mVertex* v1;
        mVertex* v2;
        mEntity* e;
        
        cout<<NbEdges<<" edges, ";
        
        for(i=0;i<NbEdges;i++)
        {
            //inFile>>GEntityId;
            fscanf(in,"%d",&GEntityId);
            if (GEntityId)
            {
                //inFile>>GEntityType>>VertexId1>>VertexId2>>EntityNbConnections>>nbPts; 
                fscanf(in,"%d %d %d %d %d",&GEntityType, &VertexId1,&VertexId2,
                       &EntityNbConnections,&nbPts);
                v1 = vertices[VertexId1-1];
                v2 = vertices[VertexId2-1];
                e = M_createE(theMesh, v1, v2, theMesh->getGEntity(GEntityId-1,GEntityType,fp));
                edges.push_back(e);
                
                for(int j=0;j<nbPts;j++)
                {
                    switch(GEntityType)
                    {
                        case 1:// inFile>>u;
                            fscanf(in,"%le",&u);
                            break;
                        case 2: //inFile>>u>>v>>patch;
                            fscanf(in,"%le %le %d",&u,&v,&patch);
                            break;
                        default: break;
                    }
                }
            }
        }  // end of reading edges
        //  assert(edges.size()==NbEdges);
        
        // *******************************
        //	Read Faces
        // *******************************
        
        std::vector<mEntity*> faces;
        mEntity* f;
        mEntity* eg[3];
        cout<<NbFaces<<" faces";
        
        for(i=0;i<NbFaces;i++)
        {
            //inFile>>GEntityId;
            fscanf(in,"%d",&GEntityId);
            if(GEntityId)
            {
                //inFile>>GEntityType>>NbEdgesOnFace;
                fscanf(in,"%d %d",&GEntityType, &NbEdgesOnFace);
                if (NbEdgesOnFace != 3){
                    fclose(in);
                    throw new mException (__LINE__,__FILE__,
                                          "wrong number of edges on a face in loadSmsFile");
                }
                //inFile>>Edge1>>Edge2>>Edge3>>nbPts;
                fscanf(in,"%d %d %d %d",&Edge1,&Edge2,&Edge3,&nbPts);
                eg[0]=edges[abs(Edge1)-1];
                eg[1]=edges[abs(Edge2)-1];
                eg[2]=edges[abs(Edge3)-1];
                //      cout<<i<<"] Edge1, Edge2, Edge3 = "<<Edge1<<","<<Edge2<<","<<Edge3<<endl;
                //      eg[0]->print();
                //      eg[1]->print();
                //      eg[2]->print();
                
                f = M_createF(theMesh, 3, eg, nullptr, theMesh->getGEntity(GEntityId-1,GEntityType,fp));
                faces.push_back(f);
                
                for(int j=0;j<nbPts;j++)
                {
                    switch(GEntityType)
                    {
                        case 0: break;
                        case 1: //inFile>>u;
                            fscanf(in,"%le",&u);
                            break;
                        case 2:// inFile>>u>>v>>patch;
                            fscanf(in,"%le %le %d",&u,&v,&patch);
                            break;
                        case 3: break;
                        default: 
                            fclose(in);
                            throw new mException (__LINE__,__FILE__,
                                                  "unknown tag in sms loadSmsFile");	      
                    }
                }
            }
        } // end of reading faces
        //  assert(faces.size()==NbFaces);
        //delete [] eg;
        
        // *******************************
        //	Read Regions
        // *******************************
        if (NbRegions==0)
            cout<<endl;
        else	
            cout<<", "<<NbRegions<<" regions ...\n";
        mEntity* fc[4];
        for(i=0;i<NbRegions;i++)
        {
            //inFile>>GEntityId;
            fscanf(in,"%d",&GEntityId);
            if (GEntityId)
            {
                //inFile>>NbFacesOnRegion;
                fscanf(in,"%d",&NbFacesOnRegion);
                if(NbFacesOnRegion == 6)
                {
                    fclose(in);
                    throw new mException (__LINE__,__FILE__,
                                          "wrong number of edges on a face in loadSmsFile");				       
                }
                //    if(NbFacesOnRegion == 4)
                {
                    //inFile>>Face1>>Face2>>Face3>>Face4>>Dummy;
                    fscanf(in,"%d %d %d %d %d",&Face1,&Face2,
                           &Face3,&Face4,&Dummy);
                    fc[0]=faces[abs(Face1)-1];
                    fc[1]=faces[abs(Face2)-1];
                    fc[2]=faces[abs(Face3)-1];
                    fc[3]=faces[abs(Face4)-1];
                    
                    M_createR(theMesh,4,fc,nullptr,theMesh->getGEntity(GEntityId-1,3,fp));
                    
                }	    
            }
        } // end of reading regions
    
    fclose(in);
    }  
} // end of namespace
