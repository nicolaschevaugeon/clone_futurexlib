/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xSubMesh.h"

#include <algorithm>
#include <boost/iterator/filter_iterator.hpp>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

#include "AOMD_OwnerManager.h"
#include "ParUtil.h"
#include "mBuildAdj.h"
#include "mEdge.h"
#include "mFace.h"
#include "mRegion.h"
#include "mVertex.h"
#include "workInProgress.h"
#include "xDebug.h"
#include "xMesh.h"
#include "xRawDataExchanger.h"


using std::cout;
using std::endl;
using std::ostream;

namespace xfem
{
using AOMD::mEdge;
using AOMD::mEntity;
using AOMD::mFace;
using AOMD::mRegion;
using AOMD::mVertex;
using AOMD::ParUtil;

void printGmshVertexElem(ostream &fs, mVertex *v, int tag1, int tag2, int &k, int formatoption)
{
   if (!formatoption)
      fs << k++ << " "
         << " 15  " << tag1 << " " << tag2 << " 1 " << v->getId() << endl;
   else
      fs << k++ << " "
         << " 15  "
         << " 3 " << tag1 << " " << tag2 << " 23 " << v->getId() << endl;
}

void printGmshFace(ostream &fs, mFace *f, int tag1, int tag2, int &k, int formatoption)
{
   //  if(recur == 2)return;

   mVertex *vx[4];
   for (int i = 0; i < f->size(0); i++)
   {
      mVertex *v = (mVertex *)f->get(0, i);
      vx[i] = v;
   }

   //{  for(int i=0;i<recur;i++) fs << " "; }
   // int nbrec = (f->isAdjacencyCreated(2))?f->size(2):0;
   if (formatoption)
   {
      int typ = (f->size(0) == 3) ? 2 : 3;
      if (typ == 2)
         fs << k++ << " " << typ << " 3  " << tag1 << " " << tag2 << " " << ParUtil::Instance()->rank() + 1 << " "
            << vx[0]->getId() << " " << vx[1]->getId() << " " << vx[2]->getId() << "\n";

      else
         fs << k++ << " " << typ << " 3 " << tag1 << " " << tag2 << " " << ParUtil::Instance()->rank() + 1 << " "
            << vx[0]->getId() << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << "\n";
   }
   else
   {
      int typ = (f->size(0) == 3) ? 2 : 3;
      if (typ == 2)
         fs << k++ << " " << typ << " " << tag1 << " " << tag2 << " " << 3 << " " << vx[0]->getId() << " " << vx[1]->getId()
            << " " << vx[2]->getId() << "\n";

      else
         fs << k++ << " " << typ << " " << tag1 << " " << tag2 << " " << 4 << " " << vx[0]->getId() << " " << vx[1]->getId()
            << " " << vx[2]->getId() << " " << vx[3]->getId() << "\n";
   }

   //{for(int i=0;i<nbrec;i++) printGmshFace(fs,(mFace*)f->get(2,i),k,recur+1);}
}
void printGmshFace(ostream &fs, mFace *f, int tag1, int tag2, int &k) { printGmshFace(fs, f, tag1, tag2, k, 0); }

void printGmshRegion(ostream &fs, mRegion *f, int tag1, int tag2, int &k, int formatoption)
{
   //  if(recur == 2)return;

   mVertex *vx[8];

   for (int i = 0; i < f->size(0); i++)
   {
      mVertex *v = (mVertex *)f->get(0, i);
      vx[i] = v;
   }

   //{for(int i=0;i<recur;i++) fs << " ";}
   // int nbrec = (f->isAdjacencyCreated(3))?f->size(3):0;
   if (formatoption)
   {
      switch (f->getType())
      {
         case mEntity::PRISM:
            fs << k++ << " 6 3 " << tag1 << " " << tag2 << " " << ParUtil::Instance()->rank() + 1 << " " << vx[0]->getId() << " "
               << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << " " << vx[4]->getId() << " "
               << vx[5]->getId() << endl;
            break;
         case mEntity::HEX:
            fs << k++ << " 5 3 " << tag1 << " " << tag2 << ParUtil::Instance()->rank() + 1 << " " << vx[0]->getId() << " "
               << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << " " << vx[4]->getId() << " "
               << vx[5]->getId() << " " << vx[6]->getId() << " " << vx[7]->getId() << endl;
            break;
         case mEntity::TET:
            fs << k++ << " 4 3 " << tag1 << " " << tag2 << " " << ParUtil::Instance()->rank() + 1 << " " << vx[0]->getId() << " "
               << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << endl;
            break;
         default:
            throw;
      }
   }
   else
   {
      switch (f->getType())
      {
         case mEntity::PRISM:
            fs << k++ << " 6 " << tag1 << " " << tag2 << " "
               << " 6 " << vx[0]->getId() << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << " "
               << vx[4]->getId() << " " << vx[5]->getId() << endl;
            break;
         case mEntity::HEX:
            fs << k++ << " 5 " << tag1 << " " << tag2 << " "
               << "8 " << vx[0]->getId() << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << " "
               << vx[4]->getId() << " " << vx[5]->getId() << " " << vx[6]->getId() << " " << vx[7]->getId() << endl;
            break;
         case mEntity::TET:
            fs << k++ << " 4 " << tag1 << " " << tag2 << " "
               << "4 " << vx[0]->getId() << " " << vx[1]->getId() << " " << vx[2]->getId() << " " << vx[3]->getId() << endl;
            break;
         default:
            throw;
      }
   }

   // for(int i=0;i<nbrec;i++)
   // printGmshRegion(fs,(mRegion*)f->get(3,i),k,recur+1);
}
void printGmshRegion(ostream &fs, mRegion *f, int tag1, int tag2, int &k) { printGmshRegion(fs, f, tag1, tag2, k, 0); }

void printGmshEdge(ostream &fs, mEdge *f, int tag1, int tag2, int &k, int formatoption)
{
   //    cout << "format option "   << formatoption  <<  (ParUtil::Instance()->rank()+1)<< endl;
   //  if(recur == 2)return;
   mVertex *vx[2];
   for (int i = 0; i < f->size(0); i++)
   {
      mVertex *v = (mVertex *)f->get(0, i);
      vx[i] = v;
   }

   if (!formatoption)
   {
      fs << k++ << " 1 " << tag1 << " " << tag2 << " "
         << "2 " << vx[0]->getId() << " " << vx[1]->getId() << endl;
   }
   else
   {
      fs << k++ << " 1 3 " << tag1 << " " << tag2 << " " << ParUtil::Instance()->rank() + 1 << " " << vx[0]->getId() << " "
         << vx[1]->getId() << endl;
   }
}
void printGmshEdge(ostream &fs, mEdge *f, int tag1, int tag2, int &k) { printGmshEdge(fs, f, tag1, tag2, k, 0); }

xSubMesh::xSubMesh(const std::string _submeshname, const xMesh &_mesh)
    : submeshname(_submeshname), mesh(_mesh), _stateuptodate(true), _boundaryuptodate(true), _paralleluptodate(true)
{
}

int xSubMesh::dim() const { return (container.size(3)) ? 3 : ((container.size(2)) ? 2 : ((container.size(1)) ? 1 : 0)); }

int xSubMesh::dim_all_process() const
{

   return dim();

}

bool xSubMesh::parallelUpToDate() const
{
   int guptodate = _paralleluptodate;

   return guptodate;
}

bool xSubMesh::stateUpToDate() const
{
   int guptodate = _stateuptodate;

   return guptodate;
}

bool xSubMesh::boundaryUpToDate() const
{
   int guptodate = _boundaryuptodate;

   return guptodate;
}

std::string xSubMesh::getName() const { return submeshname; }

const xMesh &xSubMesh::getMesh() const { return mesh; }

int xSubMesh::size(int what) const { return container.size(what); }

xIter xSubMesh::begin(int what) const { return AOMD::mLeavesIterator(container.begin(what), container.end(what)); }

xIter xSubMesh::end(int what) const { return AOMD::mLeavesIterator(container.end(what), container.end(what)); }

xtool::xRange<xIter> xSubMesh::range(int what) const { return xtool::xRange<xIter>(begin(what), end(what)); }

xIterall xSubMesh::beginall(int what) const { return container.begin(what); }
xIterall xSubMesh::endall(int what) const { return container.end(what); }
mEntity *xSubMesh::find(mEntity *e) const
{
   if (tagged.getData(*e))
      return e;
   else
      return nullptr;
}

xIterall xSubMesh::beginBnd() { return boundary.begin(); }

xIterall xSubMesh::endBnd() { return boundary.end(); }

void xSubMesh::add(mEntity *e)
{
   if (!tagged.getData(*e))
   {
      // eventually add a check to see if e is in mesh.
      container.add(e);
//       AOMD::AOMD_OwnerManager *pom = mesh.getMesh().theOwnerManager;
//       if (pom->begin(e) != pom->end(e))
//       {
//          if (ownermanager.find(e) == ownermanager.end())
//          {
//             ownermanager.insert(pom->begin(e), pom->end(e));
//          }
//       }
      tagged.setData(*e) = 1;
      _stateuptodate = false;
      _boundaryuptodate = false;
      _paralleluptodate = false;
   }
}

void xSubMesh::del(mEntity *e)
{
   if (tagged.getData(*e))
   {
      container.del(e);
//       ownermanager.erase(ownermanager.lower_bound(e), ownermanager.upper_bound(e));
      tagged.deleteData(*e);
      _stateuptodate = false;
      _boundaryuptodate = false;
      _paralleluptodate = false;
   }
}

xSubMesh::~xSubMesh() {}

void xSubMesh::modifyState(int i, int j, bool state, int with)
{
   if (!state)
   {
      std::for_each(beginall(i), endall(i), AOMD::deleteAdjFunctor(j));
   }
   else if (i > j)
   {
      std::for_each(beginall(i), endall(i), xSubMesh::xCreateDownwardFunctor(j, with, const_cast<xSubMesh &>(*this)));
   }
   else if (i < j)
   {
      std::for_each(begin(j), end(j), AOMD::createUpwardFunctor(i));
   }
}

void xSubMesh::modifyAllState()
{
   if (!stateUpToDate())
   {
      modifyState(3, 2, true);
      modifyState(3, 1, true);
      modifyState(3, 0, true);
      modifyState(2, 1, true);
      modifyState(2, 0, true);
      modifyState(1, 0, true);
      modifyState(0, 1, true);
      modifyState(0, 2, true);
      modifyState(0, 3, true);
      modifyState(1, 2, true);
      modifyState(1, 3, true);
      modifyState(2, 3, true);
      _stateuptodate = true;
   }
}

void xSubMesh::updatePartitionBoundary()
{

   modifyAllState();
   _paralleluptodate = true;
// updateBoundary();

}

xSubMesh::xCreateDownwardFunctor::xCreateDownwardFunctor(int i, int j, xSubMesh &_submesh, bool f, bool up)
    : i_dim(i), j_dim(j), submesh(_submesh), force_create(f), create_upward_too(up)
{
}

/*  void xSubMesh::xCreateDownwardFunctor::operator () (mEntity *e)
    {
    //   if((e)->isAdjacencyCreated(i_dim)) e->deleteAdjacencies(i_dim);

    for(int k=0;k<e->getNbTemplates(i_dim);k++)
    {
    mEntity *t = e->getTemplate(k,i_dim,j_dim);
    mEntity *q;
    /// entity t already exists
    if ((q = submesh.find(t)))
    {
    if (e->getClassification()->dim() <
    q->getClassification()->dim() )
    {
    q -> classify (e->getClassification());
    }
    if(!(e)->isAdjacencyCreated(i_dim))
    e->add(q);
    submesh.add(q); //useless ???
    if(create_upward_too)q->appendUnique(e);
    if(q != t)delete t; //what the fuck ???
    }
    /// add the new one to the database
    else if(force_create)
    {
    if(!(e)->isAdjacencyCreated(i_dim))
    e->add(t);
    submesh.add(t);
    if(create_upward_too)t->appendUnique(e);
    }
    /// only add existing ones
    else delete t;
    }
    };*/

void xSubMesh::xCreateDownwardFunctor::operator()(mEntity *e)
{
   //   if((e)->isAdjacencyCreated(i_dim)) e->deleteAdjacencies(i_dim);
   for (int k = 0; k < e->size(i_dim); k++)
   {
      mEntity *t = e->get(i_dim, k);
      //	cout<< "add t " << t << " " << t->getType() << endl;
      if (!submesh.find(t)) submesh.add(t);
   }
}

void xSubMesh::updateBoundary()
{
   if (!boundaryUpToDate())
   {
      updatePartitionBoundary();
      boundary.clear();

      int gdim = dim_all_process();
      ;
      int bdim = gdim - 1;
      //    cout << bdim << " bdim " << endl;
      if (bdim == -1) return;

      for (xIter it = begin(bdim); it != end(bdim); ++it)
      {
         if ((*it)->size(gdim) >= 2)
         {
            if (!(find((*it)->get(gdim, 0)) && find((*it)->get(gdim, 1)))) boundary.insert(*it);
         }
//          else
//          {
//             if ((ownermanager.find(*it) == ownermanager.end())) boundary.insert(*it);
//          }
      }


      _boundaryuptodate = true;
   }
}

bool xSubMesh::isOnBoundary(mEntity *e) const { return ((boundary.find(e) != boundary.end())); }

void xSubMesh::exportGmsh(const string &filename, int formatoption, int merging) const
{
   const bool debug = xdebug_flag;
   //-- N O D   S E C T I O N  --------------------------------------------------
   // exporting vertices
   // format : id x y z
   if (debug) cout << " hello from export_sub" << endl;
   int d = dim();
   std::stringstream outfilename;
   outfilename << filename;
   int mpisize = 1;
   int myrank = 1;

   if (mpisize != 1)
   {
      switch (merging)
      {
         case 0:
         {
            outfilename << "_" << mpisize << "_" << myrank + 1 << ".msh";
            break;
         }
         case 1:
         {
            outfilename << "merged_" << mpisize << "_proc_.msh";
            break;
         }
      }
   }

   std::ofstream f;

   ParUtil::Instance()->Barrier(59, "writting0");
   // f.open(outfilename.str().c_str(), std::ios::out|std::ios::app);
   f.open(outfilename.str().c_str(), std::ios::out);

   if (formatoption)
   {
      f << "$MeshFormat" << endl;
      f << "2 0 8" << endl;
      f << "$EndMeshFormat" << endl;
   }

   f.flush();
   cout << "P" << myrank << ": " << outfilename.str() << endl;

   if (!merging)
   {
      f << "$NOE\n";
      int NV = 0;
      for (xIter it = begin(0); it != end(0); ++it)
      {
         mVertex *v = (mVertex *)(*it);
         if (v->getType() == mEntity::VERTEX) NV++;
         // else cout << "ga"   << endl;
      }

      f << NV << "\n";
      for (xIter it = begin(0); it != end(0); ++it)
      {
         mVertex *v = (mVertex *)(*it);
         if (v->getType() == mEntity::VERTEX)
            f << v->getId() << " " << v->point()(0) << " " << v->point()(1) << " " << v->point()(2) << "\n";
      }
      f << "$ENDNOE\n";

      int NBELM = 0;
      for (int dk = 0; dk < d + 1; ++dk)
      {
         for (xIter it = begin(dk); it != end(dk); ++it)
         {
            if ((GEN_type((*it)->getClassification()) == dk) || (dk == d) || (isOnBoundary(*it))) NBELM++;
         }
      }

      f << "$ELM\n";
      f << NBELM << "\n";
      int k = 1;
      for (xIter it = begin(0); it != end(0); ++it)
      {
         mVertex *vertex = (mVertex *)(*it);
         if ((GEN_type((*it)->getClassification()) == 0) || (d == 0) || (isOnBoundary(*it)))
         {
            int tag1 = GEN_tag(vertex->getClassification());
            int tag2 = tag1;
            printGmshVertexElem(f, vertex, tag1, tag2, k, formatoption);
         }
      }

      for (xIter it = begin(1); it != end(1); ++it)
      {
         mEdge *edge = (mEdge *)(*it);
         if ((GEN_type((*it)->getClassification()) == 1) || (d == 1) || (isOnBoundary(*it)))
         {
            int tag1 = GEN_tag(edge->getClassification());
            int tag2 = tag1;
            printGmshEdge(f, edge, tag1, tag2, k, formatoption);
         }
      }

      for (xIter it = begin(2); it != end(2); ++it)
      {
         mFace *face = (mFace *)(*it);
         if ((GEN_type((*it)->getClassification()) == 2) || (d == 2) || (isOnBoundary(*it)))
         {
            int tag1 = GEN_tag(face->getClassification());
            int tag2 = tag1;
            printGmshFace(f, face, tag1, tag2, k, formatoption);
         }
      }

      for (xIter it = begin(3); it != end(3); ++it)
      {
         mRegion *region = (mRegion *)(*it);
         int tag1 = GEN_tag(region->getClassification());
         int tag2 = tag1;
         printGmshRegion(f, region, tag1, tag2, k, formatoption);
      }
      f << "$ENDELM\n";
      f.flush();
      ParUtil::Instance()->Barrier(78, "writtingend");
   }

   /*

     int nbexport = (merging)?1:ParUtil::Instance()->size();


     ParUtil::Instance()->Barrier(99, "writting00");

     for (int proc = 1 ;proc <=nbexport; ++proc){
     if (merging){
     ParUtil::Instance()->Barrier(57, "writting");
     }

     f << "$NOE\n";
     int NV = 0;
     //int NPARAM = 0;
     //int NV_GV = 0;
     //int NV_MIRROR = 0;
     for(xIter it = begin_sub(0, subset);it != end_sub(0, subset); ++it)
     {
     mVertex *v = (mVertex*)(*it);
     if(v->getType() == mEntity::VERTEX) NV++;
     }
     f << NV << "\n";
     for(xIter it = begin_sub(0, subset);it != end_sub(0, subset); ++it)
     {
     mVertex *v = (mVertex*)(*it);
     if(v->getType() == mEntity::VERTEX)
     f << v->getId() << " " << v->point()(0) << " " << v->point()(1) << " " << v->point()(2) << "\n";
     }
     f << "$ENDNOE\n";



     //-- E L M   S E C T I O N --------------------------------------------------
     // exporting elements
     int NBELM = 0;
     int d = dim_sub(subset);

     int k = 1;
     if (d>=1)
     for(xIter it = begin_sub(d, subset);it != end_sub(d, subset); ++it)  NBELM++;
     if (d>=2)
     for(xIter it = begin_sub(d, subset);it != end_sub(d, subset); ++it)  NBELM++;
     if (d>=3)
     for(xIter it = begin_sub(d, subset);it != end_sub(d, subset); ++it)  NBELM++;
     f << "$ELM\n";
     f << NBELM << "\n";

     if (d>=1)
     {
     for(xIter it = begin_sub(1, subset);it != end_sub(1, subset); ++it)
     {
     mEdge *edge = (mEdge*)(*it);
     printGmshEdge (f,edge, k, formatoption);
     }
     }

     if (d>=2)
     {
     for(xIter it = begin_sub(2, subset);it != end_sub(2, subset); ++it)
     {
     mFace *face = (mFace*)(*it);
     printGmshFace (f,face, k, formatoption);
     }
     }
     if (d>=3)
     {
     for(xIter it = begin_sub(d, subset);it != end_sub(d, subset); ++it)  NBELM++;
     //print also classified face or class at the boundary of the region
     for(xIter it = begin_sub(2, subset);it != end_sub(2, subset); ++it){
     mFace * face = (mFace*)(*it);
     if ( (F_numRegions(face)==1) ||  (F_whatInType(face)==2) ) NBELM++;
     else if (F_numRegions(face)==2){
     if (!(find_sub(F_region(face, 0 ),subset)&&find_sub(F_region(face, 1),subset))) NBELM++;
     }
     }
     for(xIter it = begin_sub(0, subset);it != end_sub(0, subset); ++it){
     if ( pGEntity geom = (*it)->getClassification())
     if (geom)
     if (GEN_type(geom))
     if (V_whatInType((mVertex *)(*it))==0)  NBELM++;
     }

     for(xIter it = begin_sub(0, subset);it != end_sub(0, subset); ++it){
     mVertex * v = (mVertex *)(*it);
     if ( pGEntity geom = (*it)->getClassification())
     if (geom)
     if (GEN_type(geom))
     if (V_whatInType((mVertex *)(*it))==0)
     {
     f << ++k << " 15 " << GEN_tag(v->getClassification()) << " "  <<  v->getId() << " 1 " <<  v->getId() << endl;
     }
     }

     for(xIter it = begin_sub(d, subset);it != end_sub(d, subset); ++it){
     mRegion *region = (mRegion*)(*it);
     printGmshRegion (f,region, k, formatoption);
     }
     //print also classified face or class at the boundary of the region
     for(xIter it = begin_sub(2, subset);it != end_sub(2, subset); ++it){
     mFace * face = (mFace*)(*it);
     if ((F_whatInType(face)==2) ||  (F_numRegions(face)==1) ) printGmshFace (f,face, k, formatoption);
     // boundary of the subregion
     else if (F_numRegions(face)==2){
     if (!(find_sub(F_region(face, 0 ),subset)&&find_sub(F_region(face, 1),subset))) printGmshFace (f,face, k, formatoption);
     }
     }

     }
     }
     f << "$ENDELM\n";
     f.flush();

     if (formatoption==1){
     ParUtil::Instance()->Barrier(58, "writting2");
     }
     }

     f.close();
   */
}

void xSubMesh::compute_bounding_box(xtensor::xPoint &min, xtensor::xPoint &max) const
{
   xgeom::xBoundingBox bb;
   for (const mEntity *pv : range(0))
   {
      xtensor::xPoint p = static_cast<const AOMD::mVertex &>(*pv).point();
      bb.inclose(p);
   }
   min = bb.min;
   max = bb.max;
}

void xUnionCreator::create(const xMesh &m, const string &name)
{
   xSubMesh &sub = m.getSubMesh(name);
   for (auto sm : all)
   {
      xSubMesh &subi = m.getSubMesh(sm);
      for (int d = 0; d < 4; d++)
      {
         for (auto e : subi.range(d)) sub.add(e);
      }
   }
   return;
}

void xIntersectionCreator::create(const xMesh &m, const string &name)
{
   xSubMesh &sub = m.getSubMesh(name);
   std::vector<string>::const_iterator its;
   for (int d = 0; d < 4; d++)
   {
      its = all.begin();
      xSubMesh &subi = m.getSubMesh(*its);
      for (xIter it = subi.begin(d); it != subi.end(d); ++it)
      {
         mEntity *e = *it;
         bool to_add = true;
         its = all.begin();
         its++;
         while (its != all.end())
         {
            xSubMesh &subj = m.getSubMesh(*its);
            if (!(subj.find(e)))
            {
               to_add = false;
               break;
            }
            its++;
         }
         if (to_add) sub.add(e);
      }
   }
   return;
}

void xFirstMinusSecondCreator::create(const xMesh &m, const string &name)
{
   xSubMesh &sub = m.getSubMesh(name);
   xSubMesh &sub1 = m.getSubMesh(first);
   xSubMesh &sub2 = m.getSubMesh(second);
   for (int d = 0; d < 4; d++)
   {
      for (xIter it = sub1.begin(d); it != sub1.end(d); ++it)
      {
         mEntity *e = *it;
         if (!sub2.find(e)) sub.add(e);
      }
   }
   sub.modifyAllState();
   sub.updateBoundary();
   return;
}

xAllCreator::xAllCreator() {}

void xAllCreator::create(const xMesh &m, const string &name)
{
   xSubMesh &sub = m.getSubMesh(name);
   for (int d = 0; d < 4; ++d)
   {
      for_each(m.begin(d), m.end(d), std::bind1st(std::mem_fun(&xSubMesh::add), &sub));
   }
}

xSubMeshCreatorFilter::xSubMeshCreatorFilter(const xEntityFilter &_filter) : filter(_filter) {}

void xSubMeshCreatorFilter::create(const xMesh &mesh, const string &name)
{
   xSubMesh &sub = mesh.getSubMesh(name);
   // xAcceptAll f;

   // cout << f(*mesh.begin(2))<< std::endl;
   // cout << filter(*mesh.begin(2))<< std::endl;
   for (int d = 0; d < 4; ++d)
   {
      boost::filter_iterator<xEntityFilter, xIter> iter_begin(filter, mesh.begin(d), mesh.end(d));
      boost::filter_iterator<xEntityFilter, xIter> iter_end(filter, mesh.end(d), mesh.end(d));
      for_each(iter_begin, iter_end, std::bind1st(std::mem_fun(&xSubMesh::add), &sub));
   }
}

xCopyCreator::xCopyCreator(const std::string &_source) : source(_source) {}

void xCopyCreator::create(const xMesh &m, const string &name)
{
   xSubMesh &sub = m.getSubMesh(name);
   xSubMesh &subsource = m.getSubMesh(source);
   for (int d = 0; d < 4; ++d)
   {
      for_each(subsource.begin(d), subsource.end(d), std::bind1st(std::mem_fun(&xSubMesh::add), &sub));
   }
   sub.modifyAllState();
   sub.updatePartitionBoundary();
   sub.updateBoundary();
}

void xAddLayerCreator::create(const xMesh &m, const string &name)
{
   xCopyCreator cop(initial);
   cop.create(m, name);
   xAddLayerModifier adl(layers, dim_growth, filter);
   adl.modify(m, name);
}

void xAddLayerModifier::modify(const xMesh &m, const string &name) const
{
   const bool debug = xdebug_flag;
   xSubMesh &sub = m.getSubMesh(name);
   if (layers == 0) return;
   // now we add a couple of layers
   // based on the layer parameter
   sub.updatePartitionBoundary();
   sub.updateBoundary();
   const int d_sub = sub.dim_all_process();
   for (int l = 1; l <= layers; l++)
   {
      if (debug) cout << "proc " << ParUtil::Instance()->rank() + 1 << " start layer " << l << endl;
      std::set<mEntity *> extra_layers;
      for (xIter it = sub.begin(dim_growth); it != sub.end(dim_growth); ++it)
      {
         mEntity *e = *it;
         for (int j = 0; j < e->size(d_sub); j++)
         {
            mEntity *ee = e->get(d_sub, j);
            if (filter(ee))
            {
               sub.add(ee);
            }
         }
      }

      sub.updatePartitionBoundary();
      sub.updateBoundary();
      if (debug) cout << "proc " << ParUtil::Instance()->rank() + 1 << " end layer " << l << endl;
   }

   if (debug) cout << "proc " << ParUtil::Instance()->rank() + 1 << " end." << endl;
   return;
}

xMesh *xSubMesh::copyInNewMesh()
{
   xMesh *mesh_result = new xMesh;
   AOMD::mMesh &mmesh_result = mesh_result->getMesh();
   // Copy each node -------------
   for (xIter it = begin(0); it != end(0); it++)
      mmesh_result.createVertex((*it)->getId(), static_cast<mVertex *>(*it)->point(), (*it)->getClassification());
   // Copy each edge ---------------
   for (xIter it = begin(1); it != end(1); it++)
      mmesh_result.createEdge((((mEdge *)*it)->get(0, 0))->getId(), (((mEdge *)*it)->get(0, 1))->getId(),
                              (*it)->getClassification());
   // Copye each face  ------------
   for (xIter it = begin(2); it != end(2); it++)
   {
      mFace *face = (mFace *)*it;
      switch (face->getType())
      {
         case mEntity::TRI:
            mmesh_result.createFaceWithVertices((face->get(0, 0))->getId(), (face->get(0, 1))->getId(),
                                                (face->get(0, 2))->getId(), (*it)->getClassification());
            break;
         case mEntity::QUAD:
            mmesh_result.createFaceWithVertices((face->get(0, 0))->getId(), (face->get(0, 1))->getId(),
                                                (face->get(0, 2))->getId(), (face->get(0, 3))->getId(),
                                                (*it)->getClassification());
            break;
         default:
            throw;
      }
   }
   // Copye each volume  ------------
   for (xIter it = begin(3); it != end(3); it++)
   {
      mEntity *e = *it;
      switch (e->getType())
      {
         case mEntity::TET:
            mmesh_result.createTetWithVertices((e->get(0, 0))->getId(), (e->get(0, 1))->getId(), (e->get(0, 2))->getId(),
                                               (e->get(0, 3))->getId(), (*it)->getClassification());
            break;
         case mEntity::HEX:
            mmesh_result.createHexWithVertices((e->get(0, 0))->getId(), (e->get(0, 1))->getId(), (e->get(0, 2))->getId(),
                                               (e->get(0, 3))->getId(), (e->get(0, 4))->getId(), (e->get(0, 5))->getId(),
                                               (e->get(0, 6))->getId(), (e->get(0, 7))->getId(), (*it)->getClassification());
            break;
         default:
            throw;
      }
   }

   classifyUnclassifiedVerices(&mmesh_result);
   // mesh_result->modifyAllState();
   mmesh_result.modifyState(3, 2, true);
   mmesh_result.modifyState(3, 1, true);
   mmesh_result.modifyState(3, 0, true);
   mmesh_result.modifyState(2, 1, true);
   mmesh_result.modifyState(2, 0, true);
   mmesh_result.modifyState(1, 0, true);
   mmesh_result.modifyState(0, 1, true);
   mmesh_result.modifyState(0, 2, true);
   mmesh_result.modifyState(0, 3, true);
   mmesh_result.modifyState(1, 2, true);
   mmesh_result.modifyState(1, 3, true);
   mmesh_result.modifyState(2, 3, true);
   return mesh_result;
}
const xMesh::partman_t &xSubMesh::getPartitionManager() const
{
   assert(xtool::workInProgress());
   return mesh.getPartitionManager();
}

void xUnifySubMeshAccrossProcess::modify(const xMesh &m, const string &subname) const
{
   using message_t = std::vector<const AOMD::mEntity *>;
   std::map<int, message_t> messages;
   xSubMesh &sub = m.getSubMesh(subname);
   const auto &partmanmesh = m.getPartitionManager();
   int dim = m.dim();
   for (int d = 0; d <= dim; d++)
   {
      for (xIter it = sub.begin(d); it != sub.end(d); ++it)
      {
         auto po = partmanmesh.getConstPartitionObject(**it);
         for (auto ro : po.getRemoteObjectsCollectionRange()) messages[ro.getProcId()].push_back(ro.getObjectAddress());
      }
   }
   std::function<void(int, const message_t &)> recvfunction = [&sub](int source, const message_t &message) {
      for (const auto &pe : message) sub.add(const_cast<mEntity *>(pe));
   };

   xtool::sendRecvRawData(messages, recvfunction, partmanmesh.getComm());
}

}  // namespace xfem
