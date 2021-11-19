/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#include "xEntityToEntity.h"

#include <iostream>

#include "GEntity.h"
#include "mAOMD.h"
#include "mAttachableDataContainer.h"
#include "xDebug.h"
#include "xDomain.h"
#include "xMesh.h"
#include "xZone.h"

using namespace std;
using namespace AOMD;
namespace xfem
{
xClassifyOn::xClassifyOn(const std::string& z) : zone_name(z), zone_id(xZone::names.getId(z))
{
   const bool debug = xdebug_flag;
   if (debug) cout << "zone_id " << zone_id << "zone_name " << zone_name << endl;
}

AOMD::mEntity* xClassifyOn::operator()(AOMD::mEntity* e)
{
   xDomain::set(*e) = zone_id;
   return e;
}

AOMD::mEntity* xUpperAdjacency::operator()(AOMD::mEntity* e)
{
   const bool debug = xdebug_flag;
   if (debug) cout << "level e " << e->getLevel() << endl;
   AOMD::mEntity* ret = e->get(e->getLevel() + 1, 0);
   if (debug) cout << "e ret   " << endl;
   if (debug) ret->print();
   return ret;
}

AOMD::mEntity* xUpperAdjacencyPartition::operator()(AOMD::mEntity* e)
{
   AOMD::mEntity* ebig = e->get(e->getLevel() + 1, 0);
   set<mEntity*> set_e;
   for (int ii = 0; ii < e->size(0); ++ii) set_e.insert(e->get(0, ii));
   xPartition partition;
   xMesh::getPartition(ebig, partition);
   for (xPartition::iterator it = partition.begin(); it != partition.end(); ++it)
   {
      mEntity* esmall = *it;
      if (esmall == ebig) return ebig;
      int sz = esmall->size(e->getLevel());
      for (int i = 0; i < sz; ++i)
      {
         set<mEntity*> set_bnd;
         mEntity* bnd = esmall->get(e->getLevel(), i);
         int szn = bnd->size(0);
         for (int j = 0; j < szn; ++j)
         {
            mEntity* nbnd = bnd->get(0, j);
            mEntity* const* pnbnd_dupl = xMesh::get_const_was_duplicated_from().getData(*nbnd);
            if (pnbnd_dupl)
               set_bnd.insert(*pnbnd_dupl);
            else
               set_bnd.insert(nullptr);
         }
         if (set_e == set_bnd) return esmall;
      }
   }
   return ebig;  // when everything else fails. The partition is probably messed up at this point
}

xCreator::xCreator() {}

AOMD::mEntity* xCreator::operator()(AOMD::mEntity* e)
{
   mEntity* const* pret = xMesh::get_const_was_created_by().getData(*e);
   return pret ? *pret : nullptr;
}

xCreatorRecursive::xCreatorRecursive() {}

AOMD::mEntity* xCreatorRecursive::operator()(AOMD::mEntity* e)
{
   mEntity* tmp1 = e;
   mEntity* tmp2 = e;
   if (!xMesh::get_const_was_created_by().getData(*e)) return nullptr;
   do
   {
      tmp1 = tmp2;
      mEntity* const* ptmp2 = xMesh::get_const_was_created_by().getData(*tmp1);
      tmp2 = ptmp2 ? *ptmp2 : nullptr;
      // cout << tmp1 << " " <<  tmp2 << std::endl;
   } while (tmp2 != nullptr);
   return tmp1;
}

xUpperCreator::xUpperCreator() {}

AOMD::mEntity* xUpperCreator::operator()(AOMD::mEntity* e)
{
   const bool debug = xdebug_flag;
   if (e == nullptr) throw;
   mEntity* const* pret = xMesh::get_const_was_created_by().getData(*e);
   if (!pret) return nullptr;
   mEntity* ret = *pret;
   if (debug) cout << "level e " << e->getLevel() << endl;
   if (debug) cout << "level ret1 " << ret->getLevel() << endl;
   if (ret->getLevel() == e->getLevel()) ret = xUpperAdjacency()(ret);
   if (debug) cout << "level ret2 " << ret->getLevel() << endl;
   if (debug) cout << "e ret   " << endl;
   if (debug) ret->print();
   return ret;
}

xPartitionCreatorRecursive::xPartitionCreatorRecursive() {}

AOMD::mEntity* xPartitionCreatorRecursive::operator()(AOMD::mEntity* e)
{
   if (e == nullptr) throw;

   mEntity* const* pret = xMesh::get_is_in_partition_of().getData(*e);
   if (pret == nullptr) return e;
   return this->operator()(*pret);
}

AOMD::mEntity* xUpperCreatorFiltered::operator()(AOMD::mEntity* e)
{
   const bool debug = xdebug_flag;
   if (debug) cout << "level e " << e->getLevel() << endl;

   mEntity* ret = xMesh::get_const_was_created_by().at(*e);
   int lev = ret->getLevel();
   if (lev == e->getLevel())
   {
      int nb = ret->size(lev + 1);
      mEntity* upperE;
      for (int j = 0; j < nb; ++j)
      {
         upperE = ret->get(lev + 1, j);
         if (filter(upperE))
         {
            ret = upperE;
            break;
         }
      }
   }

   if (debug) cout << "e ret   " << endl;
   if (debug) ret->print();
   return ret;
}

xUpperCreatorRecursive::xUpperCreatorRecursive(int _dim) : dim(_dim) {}

AOMD::mEntity* xUpperCreatorRecursive::operator()(AOMD::mEntity* e)
{
   mEntity* tmp1 = e;
   mEntity* tmp2 = e;
   if (!xMesh::get_const_was_created_by().getData(*e)) return nullptr;
   do
   {
      tmp1 = tmp2;
      mEntity* const* ptmp2 = xMesh::get_const_was_created_by().getData(*tmp1);
      tmp2 = ptmp2 ? *ptmp2 : nullptr;
   } while (tmp2);

   mEntity* ret = tmp1;
   if (ret == nullptr) return ret;
   for (int i = 0; i < dim - tmp1->getLevel(); ++i)
   {
      ret = xUpperAdjacency()(ret);
      // std::cout << ret << std::endl;
   }
   // std::cout << ret << std::endl;
   return ret;
}

xAttached::xAttached(unsigned int tag_) : tag(tag_){};

AOMD::mEntity* xAttached::operator()(AOMD::mEntity* e)
{
   mEntity* ret = e->getAttachedEntity(tag);
   return ret;
}

}  // namespace xfem
