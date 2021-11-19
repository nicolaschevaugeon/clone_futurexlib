/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xAlgorithm.h"



namespace xfem
{
// helper function for generateMpcWithL2Projection
int setDofNum(xValue<double>* v)
{
   xStateOfValue* s = v->getState();
   return ((static_cast<const xStateOfValueDof*>(s))->Numdof - 1);
}

string GetFileNameWithProcessRank(const string& base, const string& exte)
{
   return base + ".msh";
}

}  // namespace xfem
