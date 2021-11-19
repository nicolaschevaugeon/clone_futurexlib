/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

// header a redispatcher
#include <set>
#include <utility>

#include "xNonLocalInfoForKeysAndFcts.h"

using AOMD::mEntity;

namespace xfem
{
xNonLocalInfoForKeysAndFcts::xNonLocalInfoForKeysAndFcts() : it_slave_end(slave.end()), it_master_end(master.end()) {}

int xNonLocalInfoForKeysAndFcts::getAllKeysFromSlavesAndFree(xValKey::ids_size_t phys, xSpace::femKeys *keys,
                                                             xSpace::femKeys *other_master_keys) const
{
   xSpace::femKeys::iterator it = keys->begin();
   std::set<xValKey> masters;
   xValKeySub search_key;
   for (; it != keys->end();)
   {
      slave_container_const_iter_t it_find;
      search_key.setGeom(it->getGeom());
      search_key.setEnti(it->getEnti());
      if ((it_find = slave.find(search_key)) != it_slave_end)
      {
         it = keys->erase(it);
         xSpace::femKeys::const_iterator it_cur_slave = (it_find->second).begin();
         xSpace::femKeys::const_iterator it_cur_slave_end = (it_find->second).end();
         for (; it_cur_slave != it_cur_slave_end; ++it_cur_slave) masters.insert(*it_cur_slave);
      }
      else
         ++it;
   }
   int nb_free = keys->size();
   keys->reserve(keys->size() + masters.size());
   other_master_keys->reserve(masters.size());
   std::set<xValKey>::iterator itm = masters.begin();
   std::set<xValKey>::iterator itme = masters.end();
   for (; itm != itme; ++itm)
   {
      xValKey m(*itm);
      m.setPhys(phys);
      if (find(keys->begin(), keys->end(), m) != keys->end())
         other_master_keys->push_back(m);
      else
         keys->push_back(m);
   }

   return nb_free;
}

// nota :
// this return pair containing keys with the wrong phys id : i.e phys id from
// the space used to generate the NLI definition. For now it's not problematic as getCoefAndKeysForKey
// is only used in xSpacePolynomialExtendedLagrange::getFcts where only entity and geom id are adressed ....
const xNonLocalInfoForKeysAndFcts::femRelKeys &xNonLocalInfoForKeysAndFcts::getCoefAndKeysForKey(const xValKey &key) const
{
   xValKeySub search_key;
   search_key.setGeom(key.getGeom());
   search_key.setEnti(key.getEnti());
   master_container_const_iter_t it_find;
   if ((it_find = master.find(search_key)) != it_master_end)
   {
      return it_find->second;
   }
   else
   {
      cout << "Strange that you ask for a master key that doesn't exist !!!\n master key asked : " << key << endl;
      throw;
   }
}

void xNonLocalInfoForKeysAndFcts::addSlave(xValKey &key, xSpace::femKeys &keys)
{
   xValKeySub sub_key(key.getEnti(), key.getGeom());
   slave.insert(std::make_pair(sub_key, keys));
   it_slave_end = slave.end();
}
void xNonLocalInfoForKeysAndFcts::addMaster(xValKey &key, xNonLocalInfoForKeysAndFcts::femRelKeys &rel_keys)
{
   xValKeySub sub_key(key.getEnti(), key.getGeom());
   master.insert(std::make_pair(sub_key, rel_keys));
   it_master_end = master.end();
}

void xNonLocalInfoForKeysAndFcts::merge(xNonLocalInfoForKeysAndFcts *added)
{
#ifndef NDEBUG
   const bool debug = false;
   if (debug)
   {
      cout << "In merge\n";
      cout << "Size slave : " << slave.size() << endl;
   }
#endif
   // for slave merging is juste adding new entry in slave table
   slave_container_iter_t its = added->slave.begin();
   slave_container_iter_t itse = added->it_slave_end;
   for (; its != itse; ++its)
   {
      // here we don't care if it is really inserted or not as if
      // it existe it should be the same. Adding a test would verify
      // this last assertion
      slave.insert(*its);
#ifndef NDEBUG
      if (debug) cout << "add slave " << its->first.getEnti()->getId() << endl;
#endif
   }
   it_slave_end = slave.end();

#ifndef NDEBUG
   if (debug)
   {
      cout << "Size slave after : " << slave.size() << endl;
      cout << "Size master : " << master.size() << endl;
   }
#endif

   // for master merging is :
   //  - juste adding new entry in master table if new master key
   //  - adding new termes in master table entry if old master key
   //
   master_container_iter_t itm = added->master.begin();
   master_container_iter_t itme = added->it_master_end;
   master_container_iter_t it_find;
   for (; itm != itme; ++itm)
   {
      if ((it_find = master.find(itm->first)) != it_master_end)
      {
#ifndef NDEBUG
         if (debug)
         {
            cout << "fill master " << itm->first.getEnti()->getId() << endl;
            cout << "Size slave from master  : " << itm->second.size() << endl;
         }
#endif
         // loop
         xNonLocalInfoForKeysAndFcts::femRelKeys::iterator itms = itm->second.begin();
         xNonLocalInfoForKeysAndFcts::femRelKeys::iterator itmse = itm->second.end();
         for (; itms != itmse; ++itms)
         {
            (it_find->second).insert(*itms);
         }
#ifndef NDEBUG
         if (debug)
         {
            cout << "Size slave from master after : " << itm->second.size() << endl;
         }
#endif
      }
      else
      {
#ifndef NDEBUG
         if (debug) cout << "add master " << itm->first.getEnti()->getId() << endl;
#endif
         master.insert(*itm);
         it_master_end = master.end();
      }
   }

#ifndef NDEBUG
   if (debug) cout << "Size master after : " << master.size() << endl;
#endif
}
void xNonLocalInfoForKeysAndFcts::reserveMaster(int n)
{
#ifdef USE_C11_FEATURE
   master.reserve(n);
#endif
}

std::ostream &operator<<(std::ostream &ofs, const xNonLocalInfoForKeysAndFcts &ext)
{
   ofs << "xNonLocalInfoForKeysAndFcts : " << &ext << endl;
   ofs << "slave size : " << ext.slave.size() << endl;
   ofs << "slave :\n";
   xNonLocalInfoForKeysAndFcts::slave_container_const_iter_t its = ext.slave.begin();
   xNonLocalInfoForKeysAndFcts::slave_container_iter_t itse = ext.it_slave_end;
   for (; its != itse; ++its)
   {
      ofs << "slave with geom id = " << its->first.getGeom() << " entity : " << its->first.getEnti()
          << " (id = " << its->first.getEnti()->getId() << ", Level = " << its->first.getEnti()->getLevel() << ")\n";
      ofs << "depend on : \n";
      std::copy(its->second.begin(), its->second.end(), std::ostream_iterator<xValKey>(ofs, "\n"));
   }
   ofs << "master size : " << ext.master.size() << endl;
   xNonLocalInfoForKeysAndFcts::master_container_const_iter_t itm = ext.master.begin();
   xNonLocalInfoForKeysAndFcts::master_container_iter_t itme = ext.it_master_end;
   for (; itm != itme; ++itm)
   {
      ofs << "master with geom id = " << itm->first.getGeom() << " entity : " << itm->first.getEnti()
          << " (id = " << itm->first.getEnti()->getId() << ", Level = " << itm->first.getEnti()->getLevel() << ")\n";
      ofs << "master on : \n";
      xNonLocalInfoForKeysAndFcts::femRelKeys::const_iterator itms = itm->second.begin();
      xNonLocalInfoForKeysAndFcts::femRelKeys::const_iterator itmse = itm->second.end();
      for (; itms != itmse; ++itms)
      {
         ofs << itms->first << " with coef = " << itms->second << endl;
      }
   }
   return ofs;
}

xNonLocalInfoGeneratorForKeysAndFcts::xNonLocalInfoGeneratorForKeysAndFcts() : mesh(nullptr), dim(-1), generated(false) {}
xNonLocalInfoGeneratorForKeysAndFcts::xNonLocalInfoGeneratorForKeysAndFcts(xMesh &_mesh)
    : mesh(&_mesh), dim(_mesh.dim()), generated(false)
{
}
xNonLocalInfoGeneratorForKeysAndFcts::~xNonLocalInfoGeneratorForKeysAndFcts()
{
   if (generated) clearNonLocalInfoContainer();
}

void xNonLocalInfoGeneratorForKeysAndFcts::clearNonLocalInfoContainer()
{
   if (generated)
   {
      std::vector<xNonLocalInfoForKeysAndFcts *>::iterator itcend = nli_container.end();
      std::vector<xNonLocalInfoForKeysAndFcts *>::iterator itc = nli_container.begin();
      for (; itc != itcend; ++itc)
      {
         if (*itc) delete (*itc);
      }
      nli_container.clear();
      generated = false;
   }
}
void xNonLocalInfoGeneratorForKeysAndFcts::setMesh(xMesh &_mesh)
{
   if (generated) clearNonLocalInfoContainer();
   mesh = &_mesh;
}

}  // namespace xfem
