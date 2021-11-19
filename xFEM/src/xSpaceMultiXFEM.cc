/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xSpaceMultiXFEM.h"

using std::string;

namespace xfem
{
void xSpaceMultiXFEM::getKeysAndFcts(AOMD::mEntity* e, femKeys* keys, femFcts* fcts)
{
   femKeys key_base;
   femFcts fcts_base;
   base_space->getKeysAndFcts(e, &key_base, &fcts_base);
   const int resrv = key_base.size();
   keys->reserve(resrv);
   fcts->reserve(resrv);
   femKeys::iterator it = key_base.begin();
   femKeys::iterator itend = key_base.end();
   femFcts::iterator itf = fcts_base.begin();
   for (; it != itend; ++it, ++itf)
   {
      generator_approx->generateKeysAndFcts(e, *it, key_modifier, *itf, keys, fcts);
   }
   return;
}

void xSpaceMultiXFEM::getKeys(AOMD::mEntity* e, femKeys* keys)
{
   char i;
   const bool debug = false;
   string ext_add;
   femKeys key_base;
   base_space->getKeys(e, &key_base);
   keys->reserve(key_base.size());
   femKeys::iterator it = key_base.begin();
   femKeys::iterator itend = key_base.end();
   if (debug)
   {
      std::cout << "elt " << e->getId() << " in  getKeys xSpaceMultiXFEM\n";
      std::cout << "based on nodes";
      for (i = 0; i < e->size(0); ++i) std::cout << " " << e->get(0, i)->getId();
      std::cout << "\nNumber of base key " << key_base.size() << std::endl;
   }
   for (; it != itend; ++it)
   {
      generator_approx->generateKeys(*it, key_modifier, keys);
   }
   if (debug) std::cout << "Number of enriched key " << keys->size() << std::endl;
   return;
}

}  // namespace xfem
