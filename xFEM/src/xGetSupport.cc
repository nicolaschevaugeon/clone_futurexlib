/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

//stl
#include <algorithm>
// xfem
#include "xGetSupport.h"

//AOMD
#include "mEntity.h"


using namespace std;
using namespace AOMD;

namespace xfem
{

xGetSupportConformMesh::xGetSupportConformMesh(const int _dim) : dim(_dim){};
std::set < mEntity * >  xGetSupportConformMesh::operator()( mEntity *in) const
{
    std::set < mEntity * > support;
    if (in->getLevel() == dim) support.insert(in);
    else
    {
        int size_support = in->size(dim);
        for (int i = 0; i < size_support; ++i) support.insert(in->get(dim, i));
    }
    return support;
}

} // end namespace
