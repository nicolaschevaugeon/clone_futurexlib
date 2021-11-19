/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xMesh.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>

// aomd
#include "AOMD.h"
#include "mAOMD.h"
#include "mEntity.h"
#include "xMeshToolAOMD.h"

// xfem
#include "xDebug.h"
#include "xDomain.h"
#include "xElement.h"
#include "xEntityFilter.h"
#include "xMappingBuilderHolder.h"
#include "xMesh.h"
#include "xRegion.h"
#include "xRegularGrid.h"
#include "xSubMesh.h"

// xiterface/aomd
#include "xMshToAOMDReader.h"

// xtool
#include "workInProgress.h"
#include "xDataExchanger.h"
#include "xMacro.h"
#include "xPartitionManagerTools.h"

namespace xfem
{
using std::cerr;
using std::cout;
using std::endl;
using std::make_pair;
using std::ostream;
using namespace AOMD;
using xgeom::xBoundingBox;

xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_was_created_by()
{
   static datamanager_t<AOMD::mEntity*> was_created_by;
   return was_created_by;
}
const xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_const_was_created_by() { return get_was_created_by(); }

xMesh::datamanager_t<double>& xMesh::get_r_on_edge()
{
   static datamanager_t<double> r_on_edge;
   return r_on_edge;
}
const xMesh::datamanager_t<double>& xMesh::get_const_r_on_edge() { return get_r_on_edge(); }

xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_is_duplicated_in()
{
   static datamanager_t<AOMD::mEntity*> is_duplicated_in;
   return is_duplicated_in;
}
const xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_const_is_duplicated_in() { return get_is_duplicated_in(); }

xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_was_duplicated_from()
{
   static datamanager_t<AOMD::mEntity*> was_duplicated_from;
   return was_duplicated_from;
}
const xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_const_was_duplicated_from() { return get_was_duplicated_from(); }

xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_is_in_partition_of()
{
   static datamanager_t<AOMD::mEntity*> is_duplicated_in;
   return is_duplicated_in;
}
const xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_const_is_in_partition_of() { return get_is_in_partition_of(); }

void xMesh::setDataManagers()
{
   xDomain::zone();
   get_octree_level();
   get_is_hanging_on();
   get_is_hanging_by();
   get_down_group();
   get_bnd_group();
   get_was_created_by();
   get_r_on_edge();
   get_partition();
   get_is_duplicated_in();
   get_was_duplicated_from();
   get_is_in_partition_of();
}

xMesh::datamanager_t<xMesh>& xMesh::get_partition()
{
   static datamanager_t<xMesh> partition;
   return partition;
}

const xMesh::datamanager_t<xMesh>& xMesh::get_const_partition() { return get_partition(); }

xMesh::datamanager_t<int>& xMesh::xMesh::get_octree_level()
{
   static xMesh::datamanager_t<int> octree_level;
   return octree_level;
}

const xMesh::datamanager_t<int>& xMesh::xMesh::get_const_octree_level() { return get_octree_level(); }

xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_is_hanging_on()
{
   static xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*> is_hanging_on;
   return is_hanging_on;
}

const xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_const_is_hanging_on() { return get_is_hanging_on(); }

xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_is_hanging_by()
{
   static xMesh::datamanager_t<AOMD::mEntity*> is_hanging_by;
   return is_hanging_by;
}

const xMesh::datamanager_t<AOMD::mEntity*>& xMesh::get_const_is_hanging_by() { return get_is_hanging_by(); }

xMesh::datamanager_t<std::vector<AOMD::mEntity*>>& xMesh::get_down_group()
{
   static xMesh::datamanager_t<std::vector<AOMD::mEntity*>> down_group;
   return down_group;
}

const xMesh::datamanager_t<std::vector<AOMD::mEntity*>>& xMesh::get_const_down_group() { return get_down_group(); }

xMesh::datamanager_t<std::vector<AOMD::mEntity*>>& xMesh::get_bnd_group()
{
   static xMesh::datamanager_t<std::vector<AOMD::mEntity*>> bnd_group;
   return bnd_group;
}

const xMesh::datamanager_t<std::vector<AOMD::mEntity*>>& xMesh::get_const_bnd_group() { return get_bnd_group(); }

// extension de Aomd pour prendre en compte des
// sous-domaines dans le maillage

xBoundingBox compute_bounding_box(const AOMD::mEntity& e)
{
   switch (e.getType())
   {
      case mEntity::mType::VERTEX:
         return xBoundingBox{static_cast<const AOMD::mVertex&>(e).point(), static_cast<const mVertex&>(e).point()};
      default:
         xtensor::xPoint p = static_cast<const AOMD::mVertex*>(e.get(0, 0))->point();
         xgeom::xBoundingBox bb{p, p};
         for (int i = 1; i < e.size(0); ++i) bb.inclose(static_cast<const AOMD::mVertex*>(e.get(0, i))->point());
         return bb;
   }
}

xMesh::~xMesh() { clear(); }

void xMesh::clear()
{
   deleteAllSubMesh();
   for (int i = 0; i <= 3; ++i)
      for (mEntity* e : range(i)) clear(e);
   part_man.clear();
   clearOctree();
   clearGrid();
   mesh.cleanup(dim());
}

void xMesh::clear(mEntity* pe)
{
   // remove e from all the subsets
   for (auto& s_psubmesh : subsetEntities) s_psubmesh.second->del(pe);
   // remove e from all DATAMANAGER known by the xMesh
   xDomain::del(*pe);
   get_octree_level().deleteData(*pe);
   get_is_hanging_on().deleteData(*pe);
   get_is_hanging_by().deleteData(*pe);
   get_down_group().deleteData(*pe);
   get_bnd_group().deleteData(*pe);
   get_was_created_by().deleteData(*pe);
   get_r_on_edge().deleteData(*pe);
   get_partition().deleteData(*pe);
   get_is_duplicated_in().deleteData(*pe);
   get_was_duplicated_from().deleteData(*pe);
   get_is_in_partition_of().deleteData(*pe);
   // remove e from all partitition manager
   part_man.remove(*pe);
   // we should also clean up the datastructure to locate an entity (RegularGrid  and octree)
}
/* This function has been removed.
 void xMesh::del(mEntity* e)
{
   clear(e);
   mesh.del(e);
}*/

mEntity* getSource(mEntity* e)
{
   mEntity* const* pre = xMesh::get_was_duplicated_from().getData(*e);  // e->getAttachedEntity
   return pre ? getSource(*pre) : e;
}

xMesh::xMesh(MPI_Comm world, int id) : part_man(world), mesh(id) { setDataManagers(); }

xMesh::xMesh(const string& filename, MPI_Comm world) : part_man(world), mesh()
{
   setDataManagers();
   std::ifstream meshfile(filename.c_str());
   if (!meshfile.is_open())
   {
      cout << "can't open meshfile " << filename << " in xMesh Constructor " << __FILE__ << " " << __LINE__ << endl;
      throw;
   }

   xmeshtool::xMshReader(meshfile, mesh);
   meshfile.close();
   xinterface::aomd::modifyAllState(mesh);

   int proc_id;
   MPI_Comm_rank(part_man.getComm(), &proc_id);

   std::cout << "Proc " << proc_id << " nb nodes " << size(0) << " nb edges " << size(1) << " nb faces " << size(2)
             << " nb regions " << size(3) << endl;
}

xMesh::xMesh(const string& filename, xinterface::aomd::xAttachedDataManagerAOMD<int>& entities_id, MPI_Comm world)
    : part_man(world), mesh()
{
   setDataManagers();
   std::ifstream meshfile(filename.c_str());
   if (!meshfile.is_open())
   {
      cout << "can't open meshfile " << filename << " in xMesh Constructor " << __FILE__ << " " << __LINE__ << endl;
      throw;
   }
   xmeshtool::xMshReader(meshfile, mesh, entities_id);
   meshfile.close();
   xinterface::aomd::modifyAllState(mesh);
   int proc_id;
   MPI_Comm_rank(part_man.getComm(), &proc_id);

   std::cout << "Proc " << proc_id << " nb nodes " << size(0) << " nb edges " << size(1) << " nb faces " << size(2)
             << " nb regions " << size(3) << endl;
}

AOMD::mMesh& xMesh::getMesh() { return mesh; }

const AOMD::mMesh& xMesh::getMesh() const { return mesh; }

xIter xMesh::beginVertex() const { return mesh.begin(0); }

xIter xMesh::endVertex() const { return mesh.end(0); }

xIter xMesh::beginEdge() const { return mesh.begin(1); }

xIter xMesh::endEdge() const { return mesh.end(1); }

xIter xMesh::beginFace() const { return mesh.begin(2); }

xIter xMesh::endFace() const { return mesh.end(2); }

xIter xMesh::beginSolid() const { return mesh.begin(3); }

xIter xMesh::endSolid() const { return mesh.end(3); }

xIter xMesh::begin(int what) const { return mesh.begin(what); }

xIter xMesh::end(int what) const { return mesh.end(what); }

AOMD::mEntity* xMesh::find(AOMD::mEntity* pe) const { return mesh.find(pe); }

int xMesh::size(int i) const { return mesh.size(i); }

void xMesh::copyMesh(const xMesh& other)
{
   xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*> associated_new_entity;
   copyMeshInternal(other, associated_new_entity);
}

void xMesh::copyMesh(const xMesh& other, xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& associated_new_entity)
{
   copyMeshInternal(other, associated_new_entity);
}

void xMesh::copyMeshInternal(const xMesh& other,
                             xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>& associated_new_entity)
{
   // cleaning
   clear();

   // check comunicator
   MPI_Comm other_comm = other.part_man.getComm();
   MPI_Comm this_comm = part_man.getComm();
   int res;
   MPI_Comm_compare(this_comm, other_comm, &res);
   partmanAOMD_t* part_man_sub_other = nullptr;
   if (res != MPI_IDENT)
   {
      // MPI_Abort(part_man.getComm(),-123);
      // transfert local part man to new com
      part_man_sub_other = new partmanAOMD_t(this_comm);
      int dm1 = other.dim() - 1;
      xfem::xFilteredRegion<partmanAOMD_t::c_iter_object_t, xEntityFilter> fr(
          other.part_man.beginObject(), other.part_man.endObject(), [&dm1](AOMD::mEntity* f) -> bool {
             if (f->getLevel() == dm1) return true;
             return false;
          });
      xtool::createPartitionManagerForSubGroup<datamanager_t, AOMD::mEntity>(xtool::make_range(fr.begin(), fr.end()),
                                                                             other.part_man, *part_man_sub_other);
   }
   else
      part_man_sub_other = const_cast<partmanAOMD_t*>(&other.part_man);

   MPI_Comm sub_other_comm = part_man_sub_other->getComm();

   // duplication of entity and storing association in between other entity and new created entity from other entity
   for (int i = 0; i < 4; i++)
   {
      for (auto other_e : xtool::make_range(other.getMesh().beginall(i), other.getMesh().endall(i)))
      {
         mEntity* e = mesh.copyMeshEntity(other_e);
         associated_new_entity.setData(*other_e) = e;
         /* looks to be done in copyMeshEntity
         if(!find(e))
         {
          add(e);
         }
       */
      }
   }

   // set keys for exchange
   xfem::keyManagerSendAndReceive key_man(*part_man_sub_other);
   xtool::xKeyContainerSendAndRecv<keyManagerSendAndReceive::information_key_t> key_container(sub_other_comm);
   key_container.accumulateKeysAllGather(part_man_sub_other->beginObject(), part_man_sub_other->endObject(), key_man);

   // exchange/fill partition manager
   struct
   {
      typedef xtool::homogeneous_data_style_trait data_style_trait EXLIBRIS_MACRO_WARNUNUSEDTYPE;
      typedef xtool::send_and_recv_keys_communication_trait communication_trait EXLIBRIS_MACRO_WARNUNUSEDTYPE;
      typedef AOMD::mEntity* information_t;
      typedef keyManagerSendAndReceive::information_key_t information_key_t;

      information_t getInfo(information_key_t key, int sendto)
      {
         AOMD::mEntity* const* ppe = associated_new_entity->getData(*key);
         if (ppe) return *ppe;
         // there is no reason that an entity in other have not been copied to this
         else
         {
            MPI_Abort(part_man->getComm(), -1234);
            return nullptr;  // should never reach this point but covers -Wreturn-type
         }
      }
      void setInfo(information_key_t key, const information_t& info, int receivedfrom)
      {
         assert(info);

         // we used key to retrieve attached information
         AOMD::mEntity* const* ppe = associated_new_entity->getData(*key);
         assert(ppe);

         // xPartitionObject object related to this new entity
         auto po = part_man->getPartitionObject(**ppe);
         // add remote address
         po.insert(receivedfrom, info);
      }
      xMesh::partman_t* part_man;
      xinterface::aomd::xAttachedDataManagerAOMD<AOMD::mEntity*>* associated_new_entity;
   } info_man;
   info_man.part_man = &part_man;
   info_man.associated_new_entity = &associated_new_entity;
   xtool::exchangeInformation(key_container, info_man);

   if (res != MPI_IDENT) delete part_man_sub_other;

   // TODO
   if (other.regular_grid)
      createGrid();
   else
      clearGrid();
   if (other.octree_grid)
      createOctree();
   else
      clearOctree();
}
xBoundingBox xMesh::compute_bounding_box() const
{
   xBoundingBox bb;
   xtensor::xPoint& min = bb.min;
   xtensor::xPoint& max = bb.max;
   for (const mEntity* pv : range(0))
   {
      xtensor::xPoint p = static_cast<const AOMD::mVertex&>(*pv).point();
      bb.inclose(p);
   }
   // min and max must be reduced on all proc
   std::array<double, 6> mimx = {-min(0), -min(1), -min(2), max(0), max(1), max(2)};
   MPI_Allreduce(MPI_IN_PLACE, mimx.data(), 6, MPI_DOUBLE, MPI_MAX, part_man.getComm());
   min = {-mimx[0], -mimx[1], -mimx[2]};
   max = {mimx[3], mimx[4], mimx[5]};
   return xBoundingBox{min, max};
}
void xMesh::compute_bounding_box(xtensor::xPoint& min, xtensor::xPoint& max) const
{
   const auto bb = compute_bounding_box();
   min = bb.min;
   max = bb.max;
}

std::pair<double, double> xMesh::checkMesh(const int fact) const
{
   int d = dim();
   int n = 0;
   auto it = begin(d);
   auto ite = end(d);
   bool no_element = (it == ite);
   double vol, max, min, moy;
   if (no_element)
   {
      min = std::numeric_limits<double>::max();
      max = 0.;
      moy = 0.;
   }
   else
   {
      vol = xElement(*it).getVolume();
      max = vol;
      min = vol;
      moy = vol;
      for (++it, ++n; it != ite; ++it, ++n)
      {
         vol = xElement(*it).getVolume();
         if (vol > max)
            max = vol;
         else if (vol < min)
            min = vol;
         moy += vol;
      }
   }
   std::array<double, 6> mimx = {-min, max};
   MPI_Allreduce(MPI_IN_PLACE, mimx.data(), 2, MPI_DOUBLE, MPI_MAX, part_man.getComm());
   min = -mimx[0];
   max = mimx[1];
   mimx = {moy, static_cast<double>(n)};
   MPI_Allreduce(MPI_IN_PLACE, mimx.data(), 2, MPI_DOUBLE, MPI_SUM, part_man.getComm());
   assert(mimx[1] > 0.);
   moy = mimx[0] / mimx[1];
   cout << " " << endl;
   cout << "==CheckMesh output ================================================================================================"
        << endl;
   cout << "Volumes statistique : moy=" << moy << " min=" << min << " max=" << max << endl;

   for (const mEntity* pe : range(d))
   {
      const mEntity& e = *pe;
      vol = xElement(&e).getVolume();
      if (vol * fact < moy)
      {
         cout << "Element has a volume value of " << vol << " wich is " << fact << " time lesser then average value " << moy
              << endl;
         cout << "Its connectivity is :";
         for (int i = 0; i < e.size(0); ++i) cout << "  " << e.get(0, i)->getId();
         cout << endl;
      }
   }
   cout << "==End CheckMesh output ============================================================================================"
        << endl;

   return std::make_pair(min, max);
}

xSubMesh& xMesh::getSubMesh(const string& sub) const
{
   iter_subsets it = subsetEntities.find(sub);
   if (it == subsetEntities.end())
   {
      cerr << "The entities set with name " << sub << " does not exist\n";
      throw;
   }
   return *(it->second);
}

xSubMesh& xMesh::createSubMesh(const string& sub) const
{
   deleteSubMesh(sub);
   xSubMesh* subm = new xSubMesh(sub, *this);
   subsetEntities[sub] = subm;
   return *subm;
}

xSubMesh& xMesh::createSubMesh(const string& sub, xSubMeshCreator& creator) const
{
   xSubMesh* subm = &createSubMesh(sub);
   creator.create(*this, sub);
   return *subm;
}

void xMesh::deleteSubMesh(const string& sub) const
{
   iter_subsets it = subsetEntities.find(sub);
   if (it != subsetEntities.end())
   {
      delete it->second;
      subsetEntities.erase(it);
   }
}

void xMesh::deleteAllSubMesh() const
{
   for (auto& name_pmesh : subsetEntities) delete name_pmesh.second;
   subsetEntities.clear();
}

int xMesh::dim() const
{
   int dim = size(3) ? 3 : (size(2) ? 2 : (size(1) ? 1 : 0));
   assert(dim == mesh.getDim());
   return dim;
}

int xMesh::dim_global() const
{
   int dimloc = dim();
   int dimglob = dimloc;
   MPI_Allreduce(&dimloc, &dimglob, 1, MPI_INT, MPI_MAX, part_man.getComm());
   return dimglob;
}

/*
bool xMesh::isMirror(mEntity* e)
{
   mAttachableMirrorVertex* att = (mAttachableMirrorVertex*)e->getData(mirrorTag);
   if (!att) return false;
   return true;
}

*/
/*
mMirrorVertex* xMesh::getMirrorVertex(mEntity* e)
{
   mAttachableMirrorVertex* att = (mAttachableMirrorVertex*)e->getData(mirrorTag);
   if (!att) return nullptr;
   return (mMirrorVertex*)att->e;
}
*/

struct lessUpToEpsilon
{
   // given x1 and x2 it returns true
   // if x1 and x2 differ by more than epsilon
   // and x1 < x2
   // for vertices v1 and v2 it returns true
   // if the distance between v1 and v2 differ by more than epsilon
   // and x1 < x2 (alphabetique sur X, y, z)
   bool operator()(mVertex* v1, mVertex* v2) const
   {
      double eps = 1.e-06;
      return v1->point().lexicographicLessThan(v2->point(), eps);
   }
};

bool eps_equal(const double& x1, const double& x2)
{
   double eps = 1.e-6;
   if (fabs(x1 - x2) < eps) return true;
   return false;
}

/*
void xMesh::periodicAssociation()
{
   const bool debug = xdebug_flag;
   cout << "begin association" << endl;
   std::multimap<int, int> vertexAssociations;
   /// select the nodes on the boundary
   std::set<mVertex*, lessUpToEpsilon> boundary;
   typedef std::set<mVertex*, lessUpToEpsilon>::iterator bnd_iterator;
   xtensor::xPoint pmin, pmax;
   compute_bounding_box(pmin, pmax);
   std::cout << pmax << " " << pmin << std::endl;
   for (mEntity* pv : range(0))
   {
      mVertex& v = static_cast<mVertex&>(*pv);
      const double x = v.point()(0);
      const double y = v.point()(1);
      const double z = v.point()(2);
      bool onboundary;
      onboundary = (eps_equal(x, pmax(0)) || eps_equal(x, pmin(0)) || eps_equal(y, pmax(1)) || eps_equal(y, pmin(1)));
      if (getDim() == 3) onboundary = (onboundary || eps_equal(z, pmax(2)) || eps_equal(z, pmin(2)));
      if (onboundary) boundary.insert(&v);
   }
   std::cout << "size bound" << boundary.size() << std::endl;

   // group
   for (mEntity* pv : boundary)
   {
      mVertex& v = static_cast<mVertex&>(*pv);
      const double x = v.point()(0);
      const double y = v.point()(1);
      const double z = v.point()(2);

      if (eps_equal(x, pmax(0)))
      {
         mVertex vf(-1, xtensor::xPoint(pmin(0), y, z), nullptr);
         bnd_iterator itf = boundary.find(&vf);
         if (itf != boundary.end())
         {
            if (debug) cout << "the nodes " << v.getId() << " and " << (*itf)->getId() << " are opposite on x" << endl;
            if (v.getId() < (*itf)->getId())
               vertexAssociations.insert(make_pair(v.getId(), (*itf)->getId()));
            else
               vertexAssociations.insert(make_pair((*itf)->getId(), v.getId()));
         }
      }

      if (eps_equal(y, pmax(1)))
      {
         mVertex vf(-1, xtensor::xPoint(x, pmin(1), z), nullptr);
         bnd_iterator itf = boundary.find(&vf);
         if (itf != boundary.end())
         {
            if (debug) cout << "the nodes " << v.getId() << " and " << (*itf)->getId() << " are opposite on x" << endl;
            if (v.getId() < (*itf)->getId())
               vertexAssociations.insert(make_pair(v.getId(), (*itf)->getId()));
            else
               vertexAssociations.insert(make_pair((*itf)->getId(), v.getId()));
         }
      }

      if (getDim() == 3)
         if (eps_equal(z, pmax(2)))
         {
            mVertex vf(-1, xtensor::xPoint(x, y, pmin(2)), nullptr);
            bnd_iterator itf = boundary.find(&vf);
            if (itf != boundary.end())
            {
               if (debug) cout << "the nodes " << v.getId() << " and " << (*itf)->getId() << " are opposite on z" << endl;
               if (v.getId() < (*itf)->getId())
                  vertexAssociations.insert(make_pair(v.getId(), (*itf)->getId()));
               else
                  vertexAssociations.insert(make_pair((*itf)->getId(), v.getId()));
            }
         }
   }

   AOMD_Util::resolveAssociations(this, vertexAssociations);
   try
   {
      bdryLinkSetup();
   }
   catch (...)
   {
      cout << "Some problem in AOMD::mMesh::bdryLinkSetup(), catch in" << __FILE__ << ":" << __LINE__ << endl;
   }
   if (debug)
   {
      cout << "we check if evey node on the boundary has a mirror" << endl;
      cout << "if not it means we dot have a periodic mesh" << endl;
      for (mEntity* pv : boundary)
      {
         std::list<mEntity*> mirrors;
         lookupForMirrorEntities(pv, mirrors);
         cout << "mirrors.size() is " << mirrors.size() << endl;
         assert(!mirrors.empty());
      }
   }
}
*/

// if the entity e has a mirror the support is bigger!!
void xMesh::lookupSupportBasic(mEntity* pe, std::set<mEntity*>& l)
{
   const int n = pe->getLevel();
   const int d = dim();
   if (n == d)
      l.insert(pe);
   else
      for (int i = 0; i < pe->size(d); ++i) l.insert(pe->get(d, i));
   return;
}

// KO
void xMesh::lookupSupport(mEntity* pe, std::set<mEntity*>& l)
{
   lookupSupportBasic(pe, l);
   std::list<mEntity*> mirrors;
   bool domirror = false;
   if (size(3))
   {
      if (pe->getLevel() < 3)
         domirror = true;
      else if (size(2))
      {
         if (pe->getLevel() < 2) domirror = true;
      }
   }
   else if (size(1))
      if (pe->getLevel() < 1) domirror = true;

   if (domirror)
   {
      /*
            try
            {
               lookupForMirrorEntities(pe, mirrors);
            }
            catch (...)
            {
               cout << "Some problem in lookupForMirrorEntities(), catch in" << __FILE__ << ":" << __LINE__ << endl;
               cout << "e " << pe << "type  " << pe->getType() << endl;
               throw;
            }

            std::list<mEntity*>::const_iterator it;
            for (mEntity* pe : mirrors)
            {
               std::set<mEntity*> lmir;
               lookupSupportBasic(pe, lmir);
               l.insert(lmir.begin(), lmir.end());
            }
         */
   }
   return;
}

void getPartition_(mEntity* e, xPartition& partition_, xEntityFilter filter, const xMesh::datamanager_t<xMesh>& dataMesh)
{
   const bool debug = xdebug_flag;
   if (debug) cout << " In getPartition  e is ";
   if (debug) e->print();
   const xMesh* m = dataMesh.getData(*e);
   if (!m)
   {
      if (debug) cout << " pas de am " << endl;
      if (debug) cout << " filter is " << filter(e) << endl;
      if (filter(e))
      {
         partition_.insert(e);
      }
      return;
   }
   if (debug) cout << " am existe " << endl;
   for (mEntity* pe : m->range(m->dim()))
   {
      if (debug) cout << " e in partition " << endl;
      if (debug) pe->print();
      getPartition_(pe, partition_, filter, dataMesh);
   }
}

void xMesh::getPartition(mEntity* e, xPartition& partition_, xEntityFilter filter)
{
   return getPartition_(e, partition_, filter, get_const_partition());
}

void xMesh::clearOctree() const { octree_grid.reset(nullptr); }

void xMesh::createOctree() const
{
   const xgeom::xBoundingBox bb = compute_bounding_box();
   octree_grid.reset(new xgeom::xOctreeGrid(bb));
   for (auto preg : range(dim()))
   {
      std::unique_ptr<xmapping::xMapping> pmapping(xMappingBuilderHolderSingleton::instance().buildMapping(*preg));
      xgeom::xOctreeGrid::elem_data edata{preg, pmapping->eval(pmapping->COG()), pmapping->boundingBox()};
      octree_grid->addElement(edata);
   }
}

/// return the octree data structure to locate elements (xOctreeGrid).
const xgeom::xOctreeGrid& xMesh::getOctreeGrid() const
{
   if (!octree_grid) createOctree();
   return *octree_grid;
}

std::vector<std::pair<const AOMD::mEntity*, const xtensor::xPoint>> xMesh::locateElementOctree(const xtensor::xPoint& p) const
{
   if (!octree_grid) createOctree();
   std::vector<std::pair<const AOMD::mEntity*, const xtensor::xPoint>> list_pe_uvw;
   const xgeom::xOctreeGrid* pbuck = octree_grid->findBucket(p);
   if (pbuck)
   {
      for (const auto& edata : pbuck->getDatas())
      {
         if (edata.bb.contains(p))
         {
            AOMD::mEntity* pe = static_cast<AOMD::mEntity*>(edata.elem);
            std::unique_ptr<xmapping::xMapping> pmapping(xMappingBuilderHolderSingleton::instance().buildMapping(*pe));
            xtensor::xPoint uvw;
            if (pmapping->interiorCheck(p, uvw)) list_pe_uvw.push_back(std::make_pair(pe, uvw));
         }
      }
   }
   return list_pe_uvw;
}

void xMesh::createGrid() const
{
   regular_grid.reset(new xgeom::xRegularGrid(compute_bounding_box()));
   const auto& mappingbuilder = xMappingBuilderHolderSingleton::instance();
   for (const AOMD::mEntity* e : range(dim()))
   {
      std::unique_ptr<xmapping::xMapping> pmapping(mappingbuilder.buildMapping(*e));
      regular_grid->addObject(e, pmapping->boundingBox());
   }
}

void xMesh::clearGrid() const { regular_grid.reset(nullptr); }

/// return the octree data structure to locate elements (xOctreeGrid).
const xgeom::xRegularGrid& xMesh::getRegularGrid() const
{
   if (!regular_grid) createGrid();
   return *regular_grid;
}

std::vector<std::pair<const AOMD::mEntity*, const xtensor::xPoint>> xMesh::locateElement(const xtensor::xPoint& p) const
{
   if (!regular_grid) createGrid();
   std::vector<std::pair<const AOMD::mEntity*, const xtensor::xPoint>> list_e_uvw;
   const auto& mappingbuilder = xMappingBuilderHolderSingleton::instance();
   const xgeom::xRegularGrid& rg = *regular_grid;
   const xgeom::xBrick& b = rg.getBrick(p);
   for (const void* loc : b.Objects)
   {
      const AOMD::mEntity* pe = static_cast<const AOMD::mEntity*>(loc);
      std::unique_ptr<xmapping::xMapping> pmapping(mappingbuilder.buildMapping(*pe));
      xtensor::xPoint uvw;
      if (pmapping->interiorCheck(p, uvw)) list_e_uvw.push_back(std::make_pair(pe, uvw));
   }
   return list_e_uvw;
}

void xMesh::getPartitionEntity(mEntity* pe, std::set<mEntity*>& partition)
{
   const bool debug = xdebug_flag;
   if (debug) cout << " In getPartition2  e is ";
   if (debug) pe->print();
   xfem::xMesh* pm = xfem::xMesh::get_partition().getData(*pe);
   if (!pm)
   {
      if (debug) cout << " pas de am " << endl;
      partition.insert(pe);
      return;
   }
   if (pm)
   {
      if (debug) cout << " am existe " << endl;
      for (mEntity* pe : pm->range(pm->dim()))
      {
         if (debug) cout << " e in partition " << endl;
         if (debug) pe->print();
         xMesh::getPartitionEntity(pe, partition);
      }
   }
}

partmanAOMD_t& xMesh::getPartitionManager() { return part_man; }
const partmanAOMD_t& xMesh::getPartitionManager() const { return part_man; }

xSubMesh* xMesh::setBoundary()
{
   bool debug = false;
   if (debug) cout << "\nSet boundary *********\n  Dim = " << dim() << endl;
   if (boundary) deleteSubMesh("boundary");
   // boundary = new xSubMesh("boundary",*this);
   boundary = &(createSubMesh("boundary"));

   if (dim() == 0)
   {
      cout << "Warning : coded but not tested !!\n";
      boundary->add(*(this->begin(0)));
   }

   if (dim() == 1)
   {
      cout << "Warning : coded but not tested !!\n";
      for (mEntity* pv : range(0))  // pour chaque noeud
      {
         if (pv->size(1) <= 1) boundary->add(pv);
      }
   }

   if (dim() == 2)
   {
      for (mEntity* pe : range(1))  // pour chaque edge
      {
         if (debug)
         {
            cout << "  - Edge : ";
            pe->print();
         }
         if (debug) cout << "     Size upper : " << pe->size(2) << endl;
         if (pe->size(2) <= 1)
         {
            boundary->add(pe);
            if (debug) cout << "      -> Added\n";
         }
      }
   }

   if (dim() == 3)
   {
      cout << "Warning : coded but not tested !!\n";
      for (mEntity* pf : range(2))  // pour chaque face
         if (pf->size(3) <= 1) boundary->add(pf);
   }

   if (debug) cout << "End set boundary *********\n";
   boundary->modifyAllState();
   return boundary;
}

xSubMesh* xMesh::getBoundary()
{
   if (!boundary) setBoundary();
   return boundary;
}

bool xMesh::isOnBoundary(mEntity* pe)
{
   getBoundary();
   return boundary->find(pe);
}

set<mVertex*> xMesh::getNeighbors(mVertex* v)
{
   set<mVertex*> neighbors;
   for (int i = 0; i < v->size(1); i++)  // Pour chaque arete
   {
      mEntity* e = v->get(1, i);
      if (e->get(0, 0) != v)
         neighbors.insert(static_cast<mVertex*>(e->get(0, 0)));
      else
         neighbors.insert(static_cast<mVertex*>(e->get(0, 1)));
   }
   return neighbors;
}

/// cleans the partition data_manager for all entities of the current mesh.
void xMesh::cleanPartition()
{
   for (int i = dim(); i >= 0; i--)
      for (mEntity* e : range(i)) get_partition().deleteData(*e);
}

/// cleans the classification  created by xcut::xPhysSurf
void xMesh::cleanClassification()
{
   for (int i = dim(); i >= 0; i--)
      for (mEntity* e : range(i)) xDomain::del(*e);
}

}  // namespace xfem
