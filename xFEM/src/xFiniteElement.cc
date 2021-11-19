/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xFiniteElement.h"

#include <cassert>

#include "mEntity.h"
#include "xApproxFunction.h"
#include "xGeomElem.h"
#include "xSpace.h"

namespace xfem
{
using AOMD::mEntity;

xFiniteElementKeysOnly::xFiniteElementKeysOnly() : e_current(nullptr) {}

int xFiniteElementKeysOnly::sizeKey(const std::string& ityp) const { return allKeys.find(ityp)->second->size(); }

xFiniteElementKeysOnly::iterKey xFiniteElementKeysOnly ::beginKey(const std::string& ityp) { return allKeys[ityp]->begin(); }
xFiniteElementKeysOnly::iterKey xFiniteElementKeysOnly::endKey(const std::string& ityp) { return allKeys[ityp]->end(); }

xFiniteElementKeysOnly::femKeys* xFiniteElementKeysOnly::getKeys(const std::string& ityp) { return allKeys.find(ityp)->second; }
xFiniteElementKeysOnly ::~xFiniteElementKeysOnly() { clear(); }
// private function
void xFiniteElementKeysOnly::addKeys(spacePtr f, const std::string& ityp)
{
   femKeys* keys = KeysFor(ityp);
   femKeys loc;
   f->getKeys(e_current, &loc);
   keys->insert(keys->end(), loc.begin(), loc.end());
}

void xFiniteElementKeysOnly::clear()
{
   for (allKeys_t::iterator it = allKeys.begin(); it != allKeys.end(); ++it) delete it->second;
   allKeys.clear();
}

xFiniteElementKeysOnly::femKeys* xFiniteElementKeysOnly::KeysFor(const std::string& ityp)
{
   allKeys_t::iterator it = allKeys.find(ityp);
   if (it != allKeys.end()) return it->second;
   femKeys* n = new femKeys;
   allKeys.insert(make_pair(ityp, n));
   return n;
}

/////

xFiniteElement::xFiniteElement() : e_current(nullptr) {}

xFiniteElement ::~xFiniteElement() { clear(); }

void xFiniteElement::clear()
{
   for (allFcts_t::iterator it = allFcts.begin(); it != allFcts.end(); ++it) delete it->second;
   allFcts.clear();
   for (allKeys_t::iterator it = allKeys.begin(); it != allKeys.end(); ++it) delete it->second;
   allKeys.clear();
}

// void xFiniteElement::setEntity(mEntity* e)
// {
//   if (e_current == e) return;
//   clear();
//   e_current = e;
// }

int xFiniteElement::sizeKey(const std::string& ityp) const { return allKeys.find(ityp)->second->size(); }

xFiniteElement::iterKey xFiniteElement ::beginKey(const std::string& ityp) { return allKeys[ityp]->begin(); }
xFiniteElement::iterKey xFiniteElement::endKey(const std::string& ityp) { return allKeys[ityp]->end(); }

xFiniteElement::iterFct xFiniteElement ::beginFct(const std::string& ityp) { return allFcts[ityp]->begin(); }
xFiniteElement::iterFct xFiniteElement::endFct(const std::string& ityp) { return allFcts[ityp]->end(); }

xFiniteElement::femFcts* xFiniteElement::getFcts(const std::string& ityp) { return allFcts.find(ityp)->second; }
xFiniteElement::femKeys* xFiniteElement::getKeys(const std::string& ityp) { return allKeys.find(ityp)->second; }

xFiniteElement::femFcts* xFiniteElement::FctsFor(const std::string& ityp)
{
   allFcts_t::iterator it = allFcts.find(ityp);
   if (it != allFcts.end()) return it->second;
   femFcts* n = new femFcts;
   allFcts.insert(make_pair(ityp, n));
   return n;
}

xFiniteElement::femKeys* xFiniteElement::KeysFor(const std::string& ityp)
{
   allKeys_t::iterator it = allKeys.find(ityp);
   if (it != allKeys.end()) return it->second;
   femKeys* n = new femKeys;
   allKeys.insert(make_pair(ityp, n));
   return n;
}

void xFiniteElement::setKeys(mEntity* e, spacePtr f, const std::string& ityp)
{
   if (e_current != e)
   {
      clear();
      e_current = e;
   }
   femKeys* keys = KeysFor(ityp);
   keys->clear();
   addKeys(f, ityp);
}

void xFiniteElement::setKeysAndFcts(mEntity* e, spacePtr f, const std::string& ityp)
{
   if (e_current != e)
   {
      clear();
      e_current = e;
   }
   femKeys* keys = KeysFor(ityp);
   keys->clear();
   femFcts* fcts = FctsFor(ityp);
   fcts->clear();
   addKeysAndFcts(f, ityp);
}

// private function
void xFiniteElement::addKeys(spacePtr f, const std::string& ityp)
{
   femKeys* keys = KeysFor(ityp);
   femKeys loc;
   f->getKeys(e_current, &loc);
   keys->insert(keys->end(), loc.begin(), loc.end());
}

void xFiniteElement::addKeysAndFcts(spacePtr f, const std::string& ityp)
{
   femKeys* keys = KeysFor(ityp);
   femFcts* fcts = FctsFor(ityp);
   femKeys loc_keys;
   femFcts loc_fcts;
   f->getKeysAndFcts(e_current, &loc_keys, &loc_fcts);
   keys->insert(keys->end(), loc_keys.begin(), loc_keys.end());
   fcts->insert(fcts->end(), loc_fcts.begin(), loc_fcts.end());
}

}  // namespace xfem
