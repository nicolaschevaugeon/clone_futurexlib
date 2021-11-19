/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include "oLevelSetOperators.h"
#include <iostream>
#include <numeric>

using namespace std;

namespace xoctree {

oLevelSetAnalyticalInterOnActiveModifier::oLevelSetAnalyticalInterOnActiveModifier(const oLevelSetAnalytical& ls1) 
    : otherLS(ls1) 
{
}
void oLevelSetAnalyticalInterOnActiveModifier::visit(oLevelSetOnActive& ls, oKeyManager::iterator begin,  oKeyManager::iterator end)
{
    const oField &field = ls.field;
    
    for(; begin != end; ++begin)
    {
        oKey * key = *begin;
		if (key->isRegular())
		{

		        double v_point = field.getVal(key);
                field.setVal(key,std::max(v_point,otherLS.getLevelSet(key->getIJK())));
		}   

    }
}
oLevelSetAnalyticalAndShiftedActiveInterOnActiveModifier::oLevelSetAnalyticalAndShiftedActiveInterOnActiveModifier(const oLevelSetAnalytical& ls1, const oLevelSetOnActive& ls2, const double _shift) 
    : otherAnalyticalLS(ls1),otherActiveLS(ls2),shift(_shift) 
{
}
void oLevelSetAnalyticalAndShiftedActiveInterOnActiveModifier::visit(oLevelSetOnActive& ls, oKeyManager::iterator begin,  oKeyManager::iterator end)
{
    const oField &field = ls.field;
    const oField &otherActiveField = otherActiveLS.field;
    
    for(; begin != end; ++begin)
    {
        oKey * key = *begin;
		if (key->isRegular())
		{
                double v_point = otherAnalyticalLS.getLevelSet(key->getIJK());
		        double v_ls = otherActiveField.getVal(key)+shift;
                field.setVal(key,std::min(-v_point,v_ls));
		}   

    }
}


} // end namespace
