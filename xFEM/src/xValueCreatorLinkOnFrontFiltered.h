/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _VALUECREATORLINKONFRONTFILTERED_H_
#define _VALUECREATORLINKONFRONTFILTERED_H_
#include "xEntityFilter.h"
#include "xEntityToEntity.h"
#include "xValue.h"
#include "xValueCreators.h"
#include <map>
#include <vector>

namespace AOMD
{
class mVertex;
}
namespace xlinalg
{
class xCSRMatrix;
class xGraphMatrix;
}

namespace xfem
{
template <typename VT> class xValueManagerDist;
template <typename VT> class xFieldGeneric;

///// \brief Does the same as xValueCreatorLinkOnFront but filters modes that have
///// a too small support. The threshold is given with the argument "tol_" which
///// should be equal to some "h", the element characteristic length.  It has to
///// be used exactly as xValueCreatorLinkOnFront
class xValueCreatorLinkOnFrontFiltered
{
public:

    typedef std::multimap<AOMD::mVertex*, AOMD::mVertex*> Table;
    template <typename ITERFRONT>
    xValueCreatorLinkOnFrontFiltered(const ITERFRONT& begin_, const ITERFRONT& end_,
                                     xValueManagerDist<double>* v_,
                                     const double tol_, bool link_isolated_=true,
                                     xEntityFilter filt_=xAcceptAll(), xEntityToEntity entity_to_entity_=xCreator());
    ~xValueCreatorLinkOnFrontFiltered();
    xValue<double>* operator()(const xValKey& key) const;
    std::vector<AOMD::mVertex*>::iterator beginVertexIter();
    std::vector<AOMD::mVertex*>::iterator endVertexIter();
private:
    void buildTable(const int, const double, const xField<double>&, const xlinalg::xCSRMatrix&, xlinalg::xGraphMatrix&);
    xValueCreatorLinkOnFront value_creator;
    Table table;
    xValueManagerDist<double>* double_manager;
};
}

#include "xValueCreatorLinkOnFrontFiltered_imp.h"
#endif

