/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
// std
#include <numeric>
// Trellis
#include "AOMD_OwnerManager.h"
#include "ParUtil.h"
#include "mEdge.h"
// xfem
#include "xGatherMesh.h"
#include "xMesh.h"

#ifdef PARALLEL
#include "autopack.h"
#include "mpi.h"
#endif

using AOMD::mEntity;
using AOMD::mVertex;
using namespace xfem;

void BroadCastMesh(xMesh *localMesh, xMesh *globalMesh, xEntityCopyCallBack *callback)
{
   // const unsigned int  owner_tag =  AOMD::AOMD_Util::Instance()->newMeshDataId("owner");
   // const unsigned int  global_id_tag =  AOMD::AOMD_Util::Instance()->newMeshDataId("global_id");
   const unsigned int owner_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("owner");
   const unsigned int global_id_tag = AOMD::AOMD_Util::Instance()->lookupMeshDataId("global_id");
   int myrank = 0;
   int mpisize = 1;
#ifdef PARALLEL
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
#endif
   AOMD::mMesh &localmMesh = localMesh->getMesh();
   AOMD::mMesh &globalmMesh = globalMesh->getMesh();

   AOMD::AOMD_OwnerManager *pom = localmMesh.theOwnerManager;

   for (AOMD::AOMD_OwnerManager::iter it = pom->begin(); it != pom->end(); ++it)
   {
      (*it).first->attachInt(owner_tag, myrank);
   }

   for (AOMD::AOMD_OwnerManager::iter it = pom->begin(); it != pom->end(); ++it)
   {
      int cur_owner = (*it).first->getAttachedInt(owner_tag);
      int copy_pid = (*it).second.pid();
      if (cur_owner > copy_pid) (*it).first->attachInt(owner_tag, copy_pid);
   }

   std::set<mVertex *> node_on_partition_boundary;
   std::set<mVertex *> node_on_partition_boundary_owned;

   std::set<AOMD::mEdge *> edge_on_partition_boundary;
   std::set<AOMD::mEdge *> edge_on_partition_boundary_owned;

   for (AOMD::AOMD_OwnerManager::iter it = pom->begin(); it != pom->end(); ++it)
   {
      if ((*it).first->getLevel() == 0)
      {
         node_on_partition_boundary.insert(dynamic_cast<mVertex *>((*it).first));
         int owner = (*it).first->getAttachedInt(owner_tag);
         if (owner == myrank)
         {
            node_on_partition_boundary_owned.insert(dynamic_cast<mVertex *>((*it).first));
         }
      }
   }
   //  std::cout<< "size on boundary :" << 	node_on_partition_boundary_owned.size() << std::endl;

   for (AOMD::AOMD_OwnerManager::iter it = pom->begin(); it != pom->end(); ++it)
   {
      if ((*it).first->getLevel() == 1)
      {
         edge_on_partition_boundary.insert(dynamic_cast<AOMD::mEdge *>((*it).first));
         int owner = (*it).first->getAttachedInt(owner_tag);
         if (owner == myrank)
         {
            edge_on_partition_boundary_owned.insert(dynamic_cast<AOMD::mEdge *>((*it).first));
         }
      }
   }

   std::vector<int> nbnodes_global(mpisize);
   std::vector<int> nbedges_global(mpisize);

   int nbnodes_local = localmMesh.size(0) - node_on_partition_boundary.size() + node_on_partition_boundary_owned.size();
   int nbedges_local = localmMesh.size(1) - edge_on_partition_boundary.size() + edge_on_partition_boundary_owned.size();

   for (int i = 0; i < mpisize; ++i)
   {
      int nbnodes_local_i = nbnodes_local;
      int nbedges_local_i = nbedges_local;
#ifdef PARALLEL
      MPI_Bcast(&nbnodes_local_i, 1, MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(&nbedges_local_i, 1, MPI_INT, i, MPI_COMM_WORLD);
#endif
      nbnodes_global[i] = nbnodes_local_i;
      nbedges_global[i] = nbedges_local_i;
   }

   int inode = 0;
   for (int i = 0; i < myrank; ++i)
   {
      inode += nbnodes_global[i];
   }

   for (xIter it = localmMesh.begin(0); it != localmMesh.end(0); ++it)
   {
      int owner;
      if ((*it)->getAttachedInt(owner_tag, &owner))
      {
         if (owner == myrank)
         {
            (*it)->attachInt(global_id_tag, inode);
            inode++;
         }
         else
            (*it)->attachInt(global_id_tag, 0);
      }
      else
      {
         (*it)->attachInt(global_id_tag, inode);
         inode++;
      }
   }

#ifdef PARALLEL
   xExchangeMaxOnBoundary number_exchanger(global_id_tag);
   localmMesh.exchangeDataOnPartBdrys(number_exchanger, true);
#endif
   int nbnodes_total = std::accumulate(&nbnodes_global[0], &nbnodes_global[mpisize], 0);

   std::vector<double> coord(nbnodes_total * 3);
   std::vector<int> global_numbering(nbnodes_total);
   std::vector<int> global_classification_tag(nbnodes_total);
   std::vector<int> global_classification_level(nbnodes_total);

   int inodestartproc = 0;
   // broadcast nodes
   for (int i = 0; i < mpisize; ++i)
   {
      inode = inodestartproc;
      if (i == myrank)
      {
         // std::cout << "P" << myrank << std::endl;

         xIter itV = localmMesh.begin(0);
         xIter itVend = localmMesh.end(0);
         while (itV != itVend)
         {
            mVertex *v = dynamic_cast<mVertex *>(*itV);
            int owner;
            bool onboundary = v->getAttachedInt(owner_tag, &owner);
            if ((!onboundary) || (onboundary && (myrank == owner)))
            {
               Trellis_Util::mPoint p = v->point();
               global_numbering[inode] = v->getAttachedInt(global_id_tag);
               global_classification_tag[inode] = GEN_tag(v->getClassification());
               global_classification_level[inode] = GEN_type(v->getClassification());

               for (int k = 0; k < 3; ++k)
               {
                  coord[inode * 3 + k] = p(k);
               }
               inode++;
            }
            ++itV;
         }
      }
#ifdef PARALLEL
      MPI_Bcast(&coord[inodestartproc * 3], 3 * nbnodes_global[i], MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(&global_numbering[inodestartproc], nbnodes_global[i], MPI_INT, i, MPI_COMM_WORLD);

      MPI_Bcast(&global_classification_tag[inodestartproc], nbnodes_global[i], MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(&global_classification_level[inodestartproc], nbnodes_global[i], MPI_INT, i, MPI_COMM_WORLD);
#endif
      inodestartproc += nbnodes_global[i];
   }

   // Create Nodes
   std::vector<mVertex *> vertices(nbnodes_total);
   for (int i = 0; i < nbnodes_total; ++i)
   {
      int id = global_numbering[i];
      // cout << id << " " << coord[3*i] << " "  <<coord[3*i+1] << " "  << coord[3*i+2] << endl;
      int gtag = global_classification_tag[i];
      int glevel = global_classification_level[i];
      //  cout << id << " " << glevel << std::endl;
      vertices[id] =
          globalmMesh.createVertex(id, coord[3 * i], coord[3 * i + 1], coord[3 * i + 2], globalmMesh.getGEntity(gtag, glevel));
   }

   // std::cout << "P" << myrank << " done Nodes "<<  vertices.size() <<  std::endl;

   // BroadCast edges
   int nbedge_total = std::accumulate(&nbedges_global[0], &nbedges_global[mpisize], 0);
   std::vector<int> edges(nbedge_total * 2);
   global_classification_tag.clear();
   global_classification_tag.resize(nbedge_total);
   global_classification_level.clear();
   global_classification_level.resize(nbedge_total);

   // cout << "tot " << nbedge_total << endl;
   int iedgestart = 0;
   std::list<std::pair<mEntity *, int>> listLocalEdgeGlobalId;
   for (int i = 0; i < mpisize; ++i)
   {
      int iedge = iedgestart;
      if (i == myrank)
      {
         xIter itE = localmMesh.begin(1);
         xIter itEend = localmMesh.end(1);
         while (itE != itEend)
         {
            AOMD::mEdge *e = dynamic_cast<AOMD::mEdge *>((*itE));
            int owner;
            bool onboundary = e->getAttachedInt(owner_tag, &owner);
            if ((!onboundary) || (onboundary && (myrank == owner)))
            {
               listLocalEdgeGlobalId.emplace_back(std::make_pair(e, iedge));
               edges[2 * iedge] = e->get(0, 0)->getAttachedInt(global_id_tag);
               edges[2 * iedge + 1] = e->get(0, 1)->getAttachedInt(global_id_tag);
               // localedgeglobalid.insert.(e, iedge);
               global_classification_tag[iedge] = GEN_tag(e->getClassification());
               global_classification_level[iedge] = GEN_type(e->getClassification());

               ++iedge;
               //	  cout << iedge << std::endl;
            }
            ++itE;
         }
      }

#ifdef PARALLEL
      MPI_Bcast(&edges[iedgestart * 2], 2 * nbedges_global[i], MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(&global_classification_tag[iedgestart], nbedges_global[i], MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(&global_classification_level[iedgestart], nbedges_global[i], MPI_INT, i, MPI_COMM_WORLD);
#endif
      iedgestart += nbedges_global[i];
   }

   std::vector<mEntity *> globalEdgeId(nbedge_total);
   for (int i = 0; i < nbedge_total; ++i)
   {
      // std:: cout << "P" << myrank << ": " << i*5+3 << std::endl;
      int n0 = edges[i * 2];
      int n1 = edges[i * 2 + 1];
      int gtag = global_classification_tag[i];
      int glevel = global_classification_level[i];
      mEntity *edge = globalmMesh.createEdge(vertices[n0], vertices[n1], globalmMesh.getGEntity(gtag, glevel));
      globalEdgeId[i] = edge;
   }

   if (callback)
   {
      std::list<std::pair<mEntity *, int>>::iterator it = listLocalEdgeGlobalId.begin();
      std::list<std::pair<mEntity *, int>>::iterator itend = listLocalEdgeGlobalId.end();
      for (; it != itend; ++it)
      {
         mEntity *localEdge = (*it).first;
         int id = (*it).second;
         mEntity *globalEdge = globalEdgeId[id];
         (*callback)(*localEdge, *globalEdge);
      }
      for (xIter itv = localmMesh.begin(0); itv != localmMesh.end(0); ++itv)
      {
         mEntity *localVertex1 = (*itv);
         int gid = localVertex1->getAttachedInt(global_id_tag);
         mEntity *globalVertex1 = vertices[gid];
         (*callback)(*localVertex1, *globalVertex1);
      }
   }

   /*for localedgeid

     get local edge
     get id
     get global edge
     callback(localedge,  globaledge);
    */

   // Broadcast nbFaceq
   std::vector<int> nbfaces_global(mpisize);
   int nbfaces_local = localmMesh.size(2);
   for (int i = 0; i < mpisize; ++i)
   {
      int nbfaces_local_i = nbfaces_local;
#ifdef PARALLEL
      MPI_Bcast(&nbfaces_local_i, 1, MPI_INT, i, MPI_COMM_WORLD);
#endif
      nbfaces_global[i] = nbfaces_local_i;
   }

   // std::cout << "P" << myrank << " " << nbfaces_global[0] << " " << nbfaces_global[1] << std::endl;

   size_t nbfaces_total = std::accumulate(&nbfaces_global[0], &nbfaces_global[mpisize], 0);

   std::cout << "P" << myrank << " done nbfaces " << nbfaces_total << std::endl;
   // broadcast faces

   std::vector<int> faces(nbfaces_total * 5);
   global_classification_tag.resize(nbfaces_total);
   global_classification_level.resize(nbfaces_total);
   int ifa = 0;
   for (int i = 0; i < mpisize; ++i)
   {
      // std::cout << "P" << myrank << "bla " <<  i << endl;

      if (i == myrank)
      {
         // std::cout << "P" << myrank << std::endl;

         xIter itF = localmMesh.begin(2);
         xIter itFend = localmMesh.end(2);
         int ifac = ifa;
         while (itF != itFend)
         {
            AOMD::mFace *f = dynamic_cast<AOMD::mFace *>((*itF));
            int nnodes = f->size(0);
            faces[ifac * 5] = nnodes;
            int k = 1;

            global_classification_tag[ifac] = GEN_tag(f->getClassification());
            //	std::cout << 	global_classification_tag[ifac] << " ";
            global_classification_level[ifac] = GEN_type(f->getClassification());
            for (; k <= nnodes; ++k)
            {
               faces[ifac * 5 + k] = f->get(0, k - 1)->getAttachedInt(global_id_tag);
            }
            ++itF;
            ifac++;
         }
      }
#ifdef PARALLEL
      int nbface_curproc = nbfaces_global[i];
      MPI_Bcast(&faces[ifa * 5], 5 * nbface_curproc, MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(&global_classification_tag[ifa], nbface_curproc, MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(&global_classification_level[ifa], nbface_curproc, MPI_INT, i, MPI_COMM_WORLD);
#endif
      ifa += nbfaces_global[i];
   }
   // std::cout << "P" << myrank << "done B face" <<  std::endl;
   std::cout.flush();
#ifdef PARALLEL
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   // Create Faces

   for (size_t i = 0; i < nbfaces_total; ++i)
   {
      // std:: cout << "P" << myrank << ": " << i*5+3 << std::endl;
      int nnodes = faces[i * 5];
      int n0 = faces[i * 5 + 1];
      int n1 = faces[i * 5 + 2];
      int n2 = faces[i * 5 + 3];
      int gid = global_classification_tag[i];
      int glevel = global_classification_level[i];

      if (nnodes == 3)
      {
         globalmMesh.createFaceWithVertices(vertices[n0], vertices[n1], vertices[n2], globalmMesh.getGEntity(gid, glevel));
      }
      else
      {
         int n3 = faces[i * 5 + 4];

         globalmMesh.createFaceWithVertices(vertices[n0], vertices[n1], vertices[n2], vertices[n3],
                                            globalmMesh.getGEntity(gid, glevel));
      }
   }
}

xExchangeMaxOnBoundary::xExchangeMaxOnBoundary(int _data_tag) : data_tag(_data_tag) {}

void *xExchangeMaxOnBoundary::AP_alloc_and_fill_buffer(mEntity *e, AOMD::AOMD_SharedInfo &si, int tag)
{
#ifdef PARALLEL
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   int buffer_size = sizeof(mEntity *) + sizeof(int);
   void *buf = AP_alloc(si.pid(), this->tag(), buffer_size);
   mEntity **ebuf = reinterpret_cast<mEntity **>(buf);
   *(ebuf++) = si.getRemotePointer();
   int *ibuf = reinterpret_cast<int *>(ebuf);
   // std::cout << "P" << myrank << " send : "  <<    e->getAttachedInt(data_tag) << std::endl;
   *ibuf = e->getAttachedInt(data_tag);
   return buf;
#else
   return nullptr;
#endif
}

void xExchangeMaxOnBoundary::receiveData(int pid, void *buf)
{
#ifdef PARALLEL
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   mEntity **ebuf = reinterpret_cast<mEntity **>(buf);
   mEntity *e = *(ebuf);
   ebuf++;
   int *data = reinterpret_cast<int *>(ebuf);
   int number = std::max(*data, e->getAttachedInt(data_tag));
   e->attachInt(data_tag, number);
   std::cout << "P" << myrank << " received : " << *data << " " << e->getAttachedInt(data_tag) << std::endl;

#endif
}

int xExchangeMaxOnBoundary::tag() const { return 36; }
