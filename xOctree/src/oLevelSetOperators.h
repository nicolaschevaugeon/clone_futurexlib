/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

// -*- C++ -*-

#ifndef _OLEVELSET_OPERATORS_H__
#define _OLEVELSET_OPERATORS_H__

#include <cstdio>
#include <vector>

#include "oLevelSet.h"

namespace xoctree {

class oLevelSetOnActiveModifier 
{
public:
  virtual void visit(oLevelSetOnActive& , oKeyManager::iterator begin,  oKeyManager::iterator end)  = 0;
  virtual ~oLevelSetOnActiveModifier() = default;;
};

class oLevelSetAnalyticalInterOnActiveModifier : public oLevelSetOnActiveModifier
{
  public:
    oLevelSetAnalyticalInterOnActiveModifier(const oLevelSetAnalytical& );
    void visit(oLevelSetOnActive& ls, oKeyManager::iterator begin,  oKeyManager::iterator end) override;
  private: 
    const  oLevelSetAnalytical& otherLS;
};

class oLevelSetAnalyticalAndShiftedActiveInterOnActiveModifier : public oLevelSetOnActiveModifier
{
  public:
    oLevelSetAnalyticalAndShiftedActiveInterOnActiveModifier(const oLevelSetAnalytical&, const oLevelSetOnActive&, const double);
    void visit(oLevelSetOnActive& ls, oKeyManager::iterator begin,  oKeyManager::iterator end) override;
  private: 
    const  oLevelSetAnalytical& otherAnalyticalLS;
    const  oLevelSetOnActive& otherActiveLS;
    const double shift;
};


} //end namespace 

#endif





