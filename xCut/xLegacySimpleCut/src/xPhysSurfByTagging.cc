/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include "xPhysSurfByTagging.h"

#include <iostream>
#include <sstream>

#include "xAOMDEntityUtil.h"
#include "xLevelSetOperators.h"
#include "xMesh.h"
#include "xPhysSurfParameter.h"
#include "xRefCutToAOMD.h"

using AOMD::mEntity;

namespace xcut
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xPhysSurfByTagging class implementation
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////// Constructor
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xPhysSurfByTagging::xPhysSurfByTagging(xfem::xLevelSet& ls_, xPhysSurfParameter param)
    : mesh(nullptr),
      mesh_bnd(nullptr),
      region(ls_.getSupport()),
      classify_in(param.getClassifyerIn()),
      classify_out(param.getClassifyerOut()),
      ls(ls_),
      fittol(param.getFittol()),
      only_higher_dim_partition(param.getPartitionParameter()),
      promotors(param.getPromotors()),
      promotor_status(0),
      partman(region.getPartitionManager()),
      exchange_key_container(partman.getComm())
{
   // init mesh
   mesh = ls.getSupport().getMesh();

   // init exchange_key_container
   xfem::keyManagerSendAndReceive key_manager(partman);
   exchange_key_container.accumulateKeysAllGather(partman.beginObject(), partman.endObject(), key_manager);

   // setting code values
   char* res = (char*)(&code_support_cover_in);
   res[0] = TE_INCOVER;
   res[1] = SCUTEW_STAT;
   res = (char*)(&code_support_cover_out);
   res[0] = TE_OUTCOVER;
   res[1] = SCUTEW_STAT;

   // create tagging , cut and classify
   construct_(param.getFitting(), param.getKeepOldPartition(), param.getRecursive());
}
/////////////////////////////////////// End Constructor

// 2 utility function used by the destructor.

// visit all the entity of  the input mesh and apply op on each
template <class MESH, class OP>
void visit(const MESH& mesh, const OP& op)
{
   for (int i : {3, 2, 1, 0})
   {
      for (auto e : mesh.range(i)) op(e);
   }
}

// visit all the entity of  the input mesh and recursivelly on all attached xmesh on them
// and apply op on each
template <class MESH, class OP>
void rvisit(const MESH& mesh, const OP& op)
{
   for (int i : {3, 2, 1, 0})
   {
      for (auto e : mesh.range(i))
      {
         op(e);
         if (xfem::xMesh* part = xfem::xMesh::get_partition().getData(*e))
         {
            rvisit(*part, op);
         }
      }
   }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Destructor
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xPhysSurfByTagging::~xPhysSurfByTagging() { clear(); }
void xPhysSurfByTagging::clear()
{
   // if debug is set to true, we test that the tags can be cleared nicely
   const bool debug = false;
   if (!entities_tag.empty())
   {
      visit(region, [this](AOMD::mEntity* e) {
         xfem::xMesh* part = xfem::xMesh::get_partition().getData(*e);
         if (this->getTag(*e) == CUT_STAT)
         {
            if (!part)
            {
               throw xRefCutToAOMDException("Strange, partition should be here", __FILE__, __LINE__, __DATE__, __TIME__);
            }
            // first check if any of the promotor already cut the entity.
            //  if it is the case we need to cclean the partition of the entity of the partition
            if (std::any_of(promotors.begin(), promotors.end(),
                            [e](xPhysSurfByTagging* promoted) { return (promoted->getTag(*e) == CUT_STAT); }))
            {
               visit(*part, [this](AOMD::mEntity* ep) {
                  if (xfem::xMesh* partp = xfem::xMesh::get_partition().getData(*ep))
                  {
                     visit(*partp, [this](AOMD::mEntity* epp) {
                        entities_tag.deleteData(*epp);
                        for (auto promoted : promotors) promoted->unPromoteStatus(epp);
                     });
                     xfem::xMesh::get_partition().deleteData(*ep);
                  }
               });
            }
            else
            {  // No previous cut, just clean the element of the partition and the partitiion.
               visit(*part, [this](AOMD::mEntity* ep) {
                  this->entities_tag.deleteData(*ep);
                  for (auto promoted : promotors) promoted->unPromoteStatus(ep);
                  xfem::xMesh* partp = xfem::xMesh::get_partition().getData(*ep);
                  if (partp)
                  {
                     throw xRefCutToAOMDException("Strange : no partition at this level", __FILE__, __LINE__, __DATE__, __TIME__);
                  }
               });
               xfem::xMesh::get_partition().deleteData(*e);
            }
         }
      });
   }
   if (mesh_bnd)
   {
      visit(*mesh_bnd, [this](AOMD::mEntity* e) {
         if (auto part_e = xfem::xMesh::get_partition().getData(*e))
         {
            visit(*part_e, [this](AOMD::mEntity* ep) {
               for (auto promoted : promotors)
               {
                  promoted->entities_tag.deleteData(*ep);
               }
               entities_tag.deleteData(*ep);
               if (auto part_ep = xfem::xMesh::get_partition().getData(*ep))
               {
                  throw xRefCutToAOMDException("Stange : no partition expected at this level", __FILE__, __LINE__, __DATE__,
                                               __TIME__);
               }
               xfem::xMesh::get_partition().deleteData(*ep);
            });
         }
         xfem::xMesh::get_partition().deleteData(*e);
         for (auto promoted : promotors)
         {
            promoted->entities_tag.deleteData(*e);
         }
      });
   }

   if (debug)
   {
      // in debug mod we visit all the entiies that might have been tagged and remove the tag,
      // to check if we do not forget anything.
      // in release it's no needed since  entities_tag.clear() will do the job.
      std::cout << "SIZE entity_tag " << entities_tag.size() << std::endl;
      auto clear_tag = [this](AOMD::mEntity* e) { this->entities_tag.deleteData(*e); };
      rvisit(region, clear_tag);
      std::cout << "SIZE entity_tag " << entities_tag.size() << std::endl;
      if (mesh_bnd) rvisit(*mesh_bnd, clear_tag);
      std::cout << "SIZE entity_tag " << entities_tag.size() << std::endl;
      for (auto promoted : promotors)
      {
         if (auto pmeshbndother = promoted->mesh_bnd)
         {
            if (pmeshbndother) rvisit(*pmeshbndother, clear_tag);
         }
      }
      std::cout << "SIZE entity_tag " << entities_tag.size() << std::endl;
      if (!entities_tag.empty())
      {
         std::cout << "Warning entities tag not properly cleaned" << std::endl;
      }
   }
   entities_tag.clear();
   support_tag.clear();
   delete mesh_bnd;
   mesh_bnd = nullptr;
}
/////////////////////////////////////// End Destructor
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Private method
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void xPhysSurfByTagging::construct_(bool fit, bool keep_old_partition_flag, bool recursive)
{
#ifdef TIMING_MONITORING
   xtool::xDeltaTime dt;
   int idt = dt.start("level set fitting");
#endif

   if (fit)
   {
      xfem::xFitToVertices fit(fittol);
      ls.accept(fit);
   }

#ifdef TIMING_MONITORING
   dt.end(idt);
   idt = dt.start("Cut and tag (seq)");
#endif
   if (mesh_bnd)
   {
      clear();
   }

   // to simplify communicator for mesh_bnd is chosen similar to mesh. But it might be more
   // optimal to consider a smaller subcomunicator (really hard to create at this level thought....)
   mesh_bnd = new xfem::xMesh(partman.getComm());

   // generate the iso-zero, cut the mesh and classify sub entity
   // keep_old_partition_flag  indicates that we keep the old partition.
   // first : construction of the object which will do the job
   xRefCutToAOMD submesh_generator(ls, entities_tag, support_tag, mesh_bnd, classify_in, classify_out, keep_old_partition_flag,
                                   recursive, promotors, only_higher_dim_partition, xfem::xMesh::get_was_created_by(),
                                   xfem::xMesh::get_r_on_edge(), xfem::xMesh::get_partition(),
                                   xfem::xMesh::get_was_duplicated_from(), xfem::xMesh::get_is_in_partition_of());
   // second : call the method
   submesh_generator.cutMeshByRef();

#ifdef TIMING_MONITORING
   dt.end(idt);
   idt = dt.start("Tag (par)");
#endif

   //  transfer of  tag_entities and tag_support for all type of entity of the frontier of process domain
   exchangeInformation(exchange_key_container, *this);

#ifdef TIMING_MONITORING
   dt.end(idt);
#endif

   // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---
   // Here is where we have to put the transfert of tempory tags to public comon tag container
   // It contain all tag of all level-set. See multilevelset  appli for a start. Has now we deal with tag with multiple status,
   // a bit container is no more posible.
   // instead a tab of char is more natural and masking an so one
   // have to be considarate as string comparaison (string base on char tag table). See if possible ...
   // Otherwise bit comparaison with mask of mask ...
   // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---

#ifdef TIMING_MONITORING
   idt = dt.start("classify");
#endif
   // classify all entity except sub entity using tagging
   // todo :  put this in cutMeshByRef or remove classification in xPhysSurfByTagging
   int i, d = mesh->dim() + 1;
   for (i = 0; i < d; ++i)
   {
      for (auto e : region.range(i))
      {
         if (supportCoversIn(e))
            classify_in(e);
         else
            classify_out(e);
      }
   }

#ifdef TIMING_MONITORING
   dt.end(idt);
   dt.print();
#endif

   return;
}

// Promotion methode
void xPhysSurfByTagging::setPromotorStatus(mEntity* e) { promotor_status = getTag(*e); }
void xPhysSurfByTagging::promoteStatus(mEntity* e)
{
   // promoting is mandatory only if entity is not already tagged
   // if it is already tagged it means that it is a sub entity generated during
   // cutting process of this xPhysSurf
   if (getTag(*e) == NO_STAT) entities_tag.setData(*e) = promotor_status;
}

void xPhysSurfByTagging::unPromoteStatus(mEntity* e) { entities_tag.deleteData(*e); }

auto xPhysSurfByTagging::getInfo(information_key_t key, int sendto) -> information_t
{
   short data;
   char* res = (char*)(&data);
   res[0] = getTag(*key);
   res[1] = getSupportTag(*key);
   return data;
}
void xPhysSurfByTagging::setInfo(information_key_t key, const information_t& info, int receivedfrom)
{
   AOMD::mEntity* e = const_cast<AOMD::mEntity*>(key);
   char res_cur_entities;
   char res_cur_support;
   bool is_rel_out = false;
   bool is_rel_in = false;
   res_cur_entities = getTag(*e);
   res_cur_support = getSupportTag(*e);

   // here entity tag is normaly consistant betwen process : nothing to do.
   // => only a partial consistant check is done for bug tracking (only when other things have to be done).
   // At some extand the transfert and check sould be removed ...

   // if the entity is allready tagged as "support cut by iso-zero Elt wise" it will stay like that whatever it was tagged in
   // other processors there is nothing to do
   if (res_cur_support & SCUTEW_STAT)
      ;
   // if the entity is allready tagged as "support cut by iso-zero " it will stay like that if it was not tagged as "support cut
   // by iso-zero Elt wise" in  other processors
   // =>  needed to check status of remote
   else if (res_cur_support & SCUT_STAT)
   {
      const char* res = reinterpret_cast<const char*>(&(info));

      // check entity status
      if (res[0] != res_cur_entities)
      {
         std::ostringstream oss;
         oss << " Bugg : Inconsitant tagging betwen process, a node already taged as " << res_cur_entities << " is tagged "
             << res[0] << " in a other process";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // it was tagged as "cut by iso-zero Elt wise" by an other process
      if (res[1] & SCUTEW_STAT)
      {
         // "cut by iso-zero Elt wise" wins => tag
         support_tag.setData(*e) = SCUTEW_STAT;
      }
   }
   // the support have "no" or "other" tagging. Depending on what happend in other process this
   // entity may have a status wich change
   else
   {
      const char* res = reinterpret_cast<const char*>(&(info));

      // check entity status
      if (res[0] != res_cur_entities)
      {
         std::ostringstream oss;
         oss << " Bugg : Inconsitant tagging betwen process, a node already taged as " << res_cur_entities << " is tagged "
             << res[0] << " in a other process";
         throw xRefCutToAOMDException(oss.str(), __FILE__, __LINE__, __DATE__, __TIME__);
      }

      // it was tagged as "cut by iso-zero Elt wise" by an other process
      if (res[1] & SCUTEW_STAT)
      {
         // "cut by iso-zero Elt wise" wins => tag
         support_tag.setData(*e) = SCUTEW_STAT;
      }
      // it was tagged as "cut by iso-zero" by an other process
      else if (res[1] & SCUT_STAT)
      {
         // "cut by iso-zero " win against "in or out touched by iso-zero" => tag
         support_tag.setData(*e) = SCUT_STAT;
      }
      // It was tagged as "in or out touched by iso-zero" by an other process
      //    => it is potentialy in this group
      else if (res[1] & TS_RELATED)
      {
         //  check status from remote
         if (res[1] & SREL_IN_STAT)
            is_rel_in = true;
         else
            is_rel_out = true;
         //  check status from this process
         if (res_cur_support & SREL_IN_STAT)
            is_rel_in = true;
         else if (res_cur_support & SREL_OUT_STAT)
            is_rel_out = true;

         //  see if entity is having simultaniously "in touched by iso-zero" and "out touched by iso-zero"  status
         //  => if yes it have to be tagged as "cut by iso-zero"
         if (is_rel_in && is_rel_out)
            support_tag.setData(*e) = SCUT_STAT;
         else if (is_rel_in)
            support_tag.setData(*e) = SREL_IN_STAT;
         else
            support_tag.setData(*e) = SREL_OUT_STAT;
      }

   }  // end else if the entity may have a status changing

   return;
}

/////////////////////////////////////// End Private methode
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Public methode
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// update classification and cutting of the domaine (when level set attached to this xPhysSurf change)
void xPhysSurfByTagging::update(bool fit, bool keep_old_partition, bool recursive)
{
   construct_(fit, keep_old_partition, recursive);
}

// give dimension of the mesh attached to this xPhysSurf
int xPhysSurfByTagging::dim() const { return mesh->dim(); }

// get the mesh or a pointeur to it
xfem::xMesh& xPhysSurfByTagging::getMesh() { return *mesh; }
const xfem::xMesh& xPhysSurfByTagging::getMesh() const { return *mesh; }
xfem::xMesh* xPhysSurfByTagging::getMeshPtr() { return mesh; }
const xfem::xMesh* xPhysSurfByTagging::getMeshPtr() const { return mesh; }

// get the iso-zero mesh
xfem::xMesh* xPhysSurfByTagging::getMesh_bnd() { return mesh_bnd; }
const xfem::xMesh* xPhysSurfByTagging::getMesh_bnd() const { return mesh_bnd; }

// Appartenance methode
//
// support
bool xPhysSurfByTagging::supportCoversIn(mEntity* e) const
{
   unsigned short int code;
   char* res = (char*)(&code);
   res[0] = getTag(*e);
   res[1] = getSupportTag(*e);
   return (code & code_support_cover_in);
}
bool xPhysSurfByTagging::supportCoversOut(mEntity* e) const
{
   unsigned short int code;
   char* res = (char*)(&code);
   res[0] = getTag(*e);
   res[1] = getSupportTag(*e);
   return (code & code_support_cover_out);
}
bool xPhysSurfByTagging::supportBoundary(mEntity* e) const { return (getSupportTag(*e) & TS_BOUNDARY); }

bool xPhysSurfByTagging::supportCutStrictly(mEntity* e) const { return (getSupportTag(*e) & TS_CUTSTRICTLY); }

bool xPhysSurfByTagging::supportCutStrictlyEltWise(mEntity* e) const { return (getSupportTag(*e) & SCUTEW_STAT); }

// atomic
bool xPhysSurfByTagging::noTouchingIn(mEntity* e) const { return (getTag(*e) & STRICT_IN_STAT); }
bool xPhysSurfByTagging::touchingIn(mEntity* e) const { return (getTag(*e) & LOOS_IN_STAT); }
bool xPhysSurfByTagging::inIsoZero(mEntity* e) const { return (getTag(*e) & ISO_ZERO_STAT); }
bool xPhysSurfByTagging::cutStrictly(mEntity* e) const { return (getTag(*e) & CUT_STAT); }
bool xPhysSurfByTagging::touchingOut(mEntity* e) const { return (getTag(*e) & LOOS_OUT_STAT); }
bool xPhysSurfByTagging::noTouchingOut(mEntity* e) const { return (getTag(*e) & STRICT_OUT_STAT); }

// comoposed
bool xPhysSurfByTagging::strictIn(mEntity* e) const { return (getTag(*e) & TE_INDOMAIN); }
bool xPhysSurfByTagging::coversIn(mEntity* e) const { return (getTag(*e) & TE_INDOMAINCUT); }
bool xPhysSurfByTagging::strictOut(mEntity* e) const { return (getTag(*e) & TE_OUTDOMAIN); }
bool xPhysSurfByTagging::coversOut(mEntity* e) const { return (getTag(*e) & TE_OUTDOMAINCUT); }

/////////////////////////////////////// End Public methode
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xPhysSurfByTagging class implementation
///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xPhysSurfByTaggingException class implementation
//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// general exception used for all xPhysSurfByTagging throw
xPhysSurfByTaggingException::xPhysSurfByTaggingException(std::string info, std::string file, int Line, std::string date,
                                                         std::string time)
{
   std::ostringstream oss;
   oss << "In file " << file << " line " << Line << " compiled " << date << " at " << time << std::endl;
   oss << "xPhysSurfByTaggingException : " << info << std::endl;
   msg = oss.str();
}
/// general exception object : destructor
xPhysSurfByTaggingException ::~xPhysSurfByTaggingException() throw() = default;

/// mandatory what method
const char* xPhysSurfByTaggingException::what() const throw() { return this->msg.c_str(); }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xPhysSurfByTaggingException class implementation
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace xcut
