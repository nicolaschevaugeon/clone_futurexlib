/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef XFEM_XGETSUPPORT_H
#define XFEM_XGETSUPPORT_H
//stl
#include <set>


namespace AOMD {
    class mEntity;
}

namespace xfem
{

class xGetSupport
{
    public:
        virtual std::set < AOMD::mEntity * > operator()(  AOMD::mEntity *in) const = 0;
};

class xGetSupportConformMesh : public xGetSupport
{
    public:
        // Dim is the dimension of the mesh entities that will constitute the support.
        xGetSupportConformMesh(const int _dim);
        // the operator expect in->get(dim, x) and in->size(dim) to return the correct answer.
        std::set < AOMD::mEntity * > operator()(  AOMD::mEntity *in) const override;
    private:
        const int dim;
};


} // end namespace

#endif
