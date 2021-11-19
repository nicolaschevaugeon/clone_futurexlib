/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XEXTENDED_H
#define _XEXTENDED_H

#include <map>
#include <utility>
#include <vector>

#include "xMesh.h"
#include "xSpace.h"
#include "xSpacePtr.h"
#include "xValKey.h"
#include "xValueCreators.h"

// nli/NLI means in this contex non local information

namespace AOMD
{
class mEntity;
}

namespace xfem
{
class xNonLocalInfoForKeysAndFcts;
std::ostream &operator<<(std::ostream &ofs, const xNonLocalInfoForKeysAndFcts &ext);

class xNonLocalInfoForKeysAndFcts
{
  public:
   typedef std::unordered_map<xValKey, double, xHashValKey, xEqualValKey> femRelKeys;

  private:
   class xValKeySub
   {
     public:
      xValKeySub() = default;
      xValKeySub(AOMD::mEntity *e_, xValKey::ids_size_t geom_id_) : e(e_), geom_id(geom_id_) {}
      AOMD::mEntity *getEnti() { return e; }
      AOMD::mEntity *getEnti() const { return e; }
      void setEnti(AOMD::mEntity *in) { e = in; }
      xValKey::ids_size_t getGeom() { return geom_id; }
      xValKey::ids_size_t getGeom() const { return geom_id; }
      void setGeom(xValKey::ids_size_t in) { geom_id = in; }

     private:
      AOMD::mEntity *e;
      xValKey::ids_size_t geom_id;
      friend int compkey(const xValKeySub &c1, const xValKeySub &c2)
      {
         if (c1.e > c2.e) return 1;
         if (c1.e < c2.e) return -1;
         if (c1.geom_id > c2.geom_id) return 1;
         if (c1.geom_id < c2.geom_id) return -1;
         return 0;
      }
   };
   struct xHashValKeySub
   {
      int operator()(const xValKeySub &key) const
      {
         AOMD::mEntity *e = key.getEnti();
         return e->getId();
      }
   };
   struct xEqualValKeySub
   {
      bool operator()(const xValKeySub &key1, const xValKeySub &key2) const { return compkey(key1, key2) == 0; }
   };

   typedef std::unordered_map<xValKeySub, femRelKeys, xHashValKeySub, xEqualValKeySub> master_container_t;
   typedef std::unordered_map<xValKeySub, xSpace::femKeys, xHashValKeySub, xEqualValKeySub> slave_container_t;

   typedef master_container_t::iterator master_container_iter_t;
   typedef slave_container_t::iterator slave_container_iter_t;
   typedef master_container_t::const_iterator master_container_const_iter_t;
   typedef slave_container_t::const_iterator slave_container_const_iter_t;
   master_container_t master;
   slave_container_t slave;
   slave_container_iter_t it_slave_end;
   master_container_iter_t it_master_end;

  public:
   xNonLocalInfoForKeysAndFcts();
   int getAllKeysFromSlavesAndFree(xValKey::ids_size_t phys, xSpace::femKeys *keys, xSpace::femKeys *other_master_keys) const;
   const femRelKeys &getCoefAndKeysForKey(const xValKey &key) const;
   void addSlave(xValKey &key, xSpace::femKeys &keys);
   void addMaster(xValKey &key, femRelKeys &keys);
   void merge(xNonLocalInfoForKeysAndFcts *added);
   void reserveMaster(int n);
   friend std::ostream &operator<<(std::ostream &ofs, const xNonLocalInfoForKeysAndFcts &ext);
};

class xNonLocalInfoGeneratorForKeysAndFcts
{
  public:
   xNonLocalInfoGeneratorForKeysAndFcts();
   xNonLocalInfoGeneratorForKeysAndFcts(xMesh &_mesh);
   virtual void generateNonLocalInfoForKeysAndFcts(spacePtr_t space, const int polynomial_order) = 0;
   virtual void clearNonLocalInfoForKeysAndFcts() = 0;
   virtual void setMesh(xMesh &mesh);
   virtual ~xNonLocalInfoGeneratorForKeysAndFcts();
   virtual const xNonLocalInfoForKeysAndFcts *getNonLocalInfo(const AOMD::mEntity &e) const = 0;

  protected:
   xMesh *mesh;
   int dim;
   bool generated;
   std::vector<xNonLocalInfoForKeysAndFcts *> nli_container;
   void clearNonLocalInfoContainer();
};

template <typename S, template <class> class DATAMANAGER = xfem::xMesh::datamanager_t>
class xNonLocalInfoGeneratorForKeysAndFctsForHangingNode : public xNonLocalInfoGeneratorForKeysAndFcts
{
  public:
   xNonLocalInfoGeneratorForKeysAndFctsForHangingNode(
       const DATAMANAGER<AOMD::mEntity *> &isHangingOn = xMesh::get_const_is_hanging_on(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup = xMesh::get_const_down_group(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup = xMesh::get_const_bnd_group(), double thresold = 1.e-12);
   xNonLocalInfoGeneratorForKeysAndFctsForHangingNode(
       xMesh &mesh, const DATAMANAGER<AOMD::mEntity *> &isHangingOn = xMesh::get_const_is_hanging_on(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup = xMesh::get_const_down_group(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup = xMesh::get_const_bnd_group(), double thresold = 1.e-12);
   xNonLocalInfoGeneratorForKeysAndFctsForHangingNode(
       xMesh &mesh, DATAMANAGER<xNonLocalInfoForKeysAndFcts *> &_nli_data,
       const DATAMANAGER<AOMD::mEntity *> &isHangingOn = xMesh::get_const_is_hanging_on(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup = xMesh::get_const_down_group(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup = xMesh::get_const_bnd_group(), double thresold = 1.e-12);
   ~xNonLocalInfoGeneratorForKeysAndFctsForHangingNode();
   void generateNonLocalInfoForKeysAndFcts(spacePtr_t space, const int polynomial_order) override;
   void clearNonLocalInfoForKeysAndFcts() override;
   const xNonLocalInfoForKeysAndFcts *getNonLocalInfo(const AOMD::mEntity &e) const override;

  private:
   const DATAMANAGER<AOMD::mEntity *> &isHangingOn;
   const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup;
   const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup;
   std::unique_ptr<DATAMANAGER<xNonLocalInfoForKeysAndFcts *>> pnli_data;
   DATAMANAGER<xNonLocalInfoForKeysAndFcts *> &nli_data;
   const double thresold;
   xValueCreator<xValueDouble> creator_reg;
   S l2_solver;
};

template <typename S, template <class> class DATAMANAGER = xfem::xMesh::datamanager_t>
class xNonLocalInfoGeneratorForKeysAndFctsForHangingNodeOldAndCurrent : public xNonLocalInfoGeneratorForKeysAndFcts
{
  public:
   xNonLocalInfoGeneratorForKeysAndFctsForHangingNodeOldAndCurrent(
       xMesh &mesh, const DATAMANAGER<AOMD::mEntity *> &isHangingOn_ = xMesh::get_const_is_hanging_on(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup_ = xMesh::get_const_down_group(),
       const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup_ = xMesh::get_const_bnd_group(), double threshold_ = 1.e-12)
       : isHangingOn(isHangingOn_),
         downGroup(downGroup_),
         bndGroup(bndGroup_),
         threshold(threshold_),
         old(nullptr),
         current(nullptr)
   {
      current =
          new xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S>(mesh, nli_data, isHangingOn, downGroup, bndGroup, threshold);
   }
   ~xNonLocalInfoGeneratorForKeysAndFctsForHangingNodeOldAndCurrent() { clearNonLocalInfoForKeysAndFcts(); }
   void setMesh(xMesh &mesh) override
   {
      if (old)
      {
         delete old;
         old = NULL;
      }
      old = current;
      current =
          new xNonLocalInfoGeneratorForKeysAndFctsForHangingNode<S>(mesh, nli_data, isHangingOn, downGroup, bndGroup, threshold);
   }
   const xNonLocalInfoForKeysAndFcts *getNonLocalInfo(const AOMD::mEntity &e) const override { return nli_data.getData(e); }
   void generateNonLocalInfoForKeysAndFcts(spacePtr_t space, const int polynomial_order) override
   {
      return current->generateNonLocalInfoForKeysAndFcts(space, polynomial_order);
   }
   void clearNonLocalInfoForKeysAndFcts() override
   {
      if (old)
      {
         delete old;
         old = NULL;
      }
      if (current)
      {
         delete current;
         current = NULL;
      }
   }

  private:
   double threshold;
   const DATAMANAGER<AOMD::mEntity *> &isHangingOn;
   const DATAMANAGER<std::vector<AOMD::mEntity *>> &downGroup;
   const DATAMANAGER<std::vector<AOMD::mEntity *>> &bndGroup;
   xNonLocalInfoGeneratorForKeysAndFcts *old;
   xNonLocalInfoGeneratorForKeysAndFcts *current;
   DATAMANAGER<xNonLocalInfoForKeysAndFcts *> nli_data;
};

}  // namespace xfem

// implementation
#include "xNonLocalInfoForKeysAndFcts_imp.h"

#endif
