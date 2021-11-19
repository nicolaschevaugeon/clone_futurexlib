/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

// -*- C++ -*-

#ifndef _OLEVELSET_H__
#define _OLEVELSET_H__

#include <cstdio>
#include <iostream>
#include <vector>

#include "oOctree.h"
#include "oKeyManager.h"
#include "oField.h"

namespace xoctree {

  class oLevelSetAnalyticalInterOnActiveModifier;

  class oLevelSet 
  {
  public:
    virtual ~oLevelSet() = default;;
    oLevelSet( const oMapping& m, int lmax);
    virtual void   getLevelSetValues (const int * ijk, int level, 
				      std::vector<double>& ls) const = 0;
    virtual double getLevelSet(const  int * ijk_node_fine) const;
    virtual void getLevelSetCurvatures(const  int * ijk, int level, std::vector<double>& curvatures) const;
    virtual double getLevelSetCurvature(const  int * ijk_node, int level) const;
    typedef std::function<double (const int*)> f_ijk_t;
    f_ijk_t getFctIJK() { return bind1st(std::mem_fun(&oLevelSet::getLevelSet), this); }
  protected:
    const oMapping& mapping;
    const int level_max;
    const int dim;
    const oTopo& topo;
  };
  
  class oLevelSetDiscrete : public  oLevelSet
  {
  public:
    oLevelSetDiscrete(const oMapping& m, int l, double* ls_);
    void   getLevelSetValues (const int * ijk, int level, std::vector<double>& ls) const override;
    double getLevelSet(const  int * ijk_node_fine) const override;
//    virtual void getLevelSetCurvatures(const  int * ijk, int level, std::vector<double>& curvatures) const ;
    double getLevelSetCurvature(const  int * ijk_node_fine, int level) const override ;///Attention: false values for nodes on the boundary: F.D. operators should be improved
    void   getLevelSetIndices (const int * ijk, int level, std::vector<int>& ls_idx) const;
    int getLevelSetIndex(const  int * ijk_node_fine) const;
  private:
    const double* ls_ijk_fine;
    const int powp1;

    /// Is the cell on the boundary ?
    /// \param[in] b : number of the bnd (1 to 4 for 2D, 1 to 6 for 3D)
    /// \param[in] ijk : ijk of the cell
    /// \param[in] level : level of the cell
    bool is_node_on_boundary(int b, const int* ijk, int level) const;//Copied from oOctree and adapted for nodes
  };

  class oLevelSetAnalytical : public  oLevelSet
  {
  public:
    typedef std::function<double (const double&, const double&, const double&)> ls_xyz_t;
    oLevelSetAnalytical(const oMapping& m, int l, ls_xyz_t ls_);
    void   getLevelSetValues (const int * ijk, int level, std::vector<double>& ls) const override;
    double getLevelSet(const  int * ijk_node_fine) const override;
//    virtual void getLevelSetCurvatures(const  int * ijk, int level, std::vector<double>& curvatures) const;
    double getLevelSetCurvature(const  int * ijk_node, int level) const override;
    void convertToDiscrete(double* ls_) const;
  private:
    ls_xyz_t ls_xyz;
  };

 class oLevelSetOnActive : public  oLevelSet
 {
   
 public:
   oLevelSetOnActive(const oOctree& octree,  const oField& f);
   
   void   getLevelSetValues (const int * ijk, int level, std::vector<double>& ls) const override;
   
   double getLevelSet(const  int * ijk_node_fine) const override;
   
   //virtual void   setLevelSet(const  int * ijk_node_fine, const double& v);
   //void refine(oOctree::const_iterator cell, int level, const int* ijk);
   //void derefine(oOctree::const_iterator cell, int level, const int* ijk);
   //void createStencils();
   

 private:
   friend class oLevelSetAnalyticalInterOnActiveModifier;
   friend class oLevelSetAnalyticalAndShiftedActiveInterOnActiveModifier;
   const oOctree& octree;
   const oField& field;
   const oKeyManager& key_manager;

   //oKeyManager& keys;
   
   //void createKeys();
   //void getKeysOnElt (const int* ijk_ori,  int level, std::vector<oKey*> & elt_keys) const;
   
   //int NbActiveNodes, NbFreeNodes;
   
 };



  class oFilterCriteriaOnSign : public oFilterCriteria
  {
    public:
      oFilterCriteriaOnSign(const oMapping& m, const oLevelSet& ls_, bool revert=false);
      bool operator()(oOctree::const_iterator cell, int level, 
                    const int * ijk, 
                    oOctree::const_iterator children_beg, 
                    oOctree::const_iterator children_end) const override;

    private:
      const oMapping& mapping;
      const oLevelSet& ls;
      const int dim;
      bool revert;
  };

  class oFilterCriteriaCrossed : public oFilterCriteria
  {
    public:
      oFilterCriteriaCrossed(const oMapping& m, const oLevelSet& ls_, bool revert=false);
      bool operator()(oOctree::const_iterator cell, int level,
                    const int * ijk,
                    oOctree::const_iterator children_beg,
                    oOctree::const_iterator children_end) const override;

    private:
      const oMapping& mapping;
      const oLevelSet& ls;
      const int dim;
      bool revert;
  };

  class oFilterCriteriaOnVolumeFraction : public oFilterCriteria
  {
    public:
      oFilterCriteriaOnVolumeFraction(const oMapping& m, const oLevelSet& ls_);
      bool operator()(oOctree::const_iterator cell, int level, 
                    const int * ijk, 
                    oOctree::const_iterator children_beg, 
                    oOctree::const_iterator children_end) const override;

    private:
      const oMapping& mapping;
      const oLevelSet& ls;
      const int dim;
  };
  
class oRefinementCriteriaOnSign : public oRefinementCriteria
{
public:
  oRefinementCriteriaOnSign(const oMapping& m, const oLevelSet& ls_);
  bool operator()(oOctree::const_iterator cell, int level, 
		  const int * ijk, 
		  oOctree::const_iterator children_beg, 
		  oOctree::const_iterator children_end) const override;

private:
  const oMapping& mapping;
  const oLevelSet& ls;
  const int dim;
};

/// Refine only cells crossed by the iso-zero, but do not derefine too much away from the interface (minimum level = levelMin)
class oRefinementCriteriaOnSignConstrained : public oRefinementCriteria
{
public:
  oRefinementCriteriaOnSignConstrained(const oMapping& m, const oLevelSet& ls_, int levelMin_, int targetLevelMax_ = 100);
  bool operator()(oOctree::const_iterator cell, int level, 
                  const int * ijk, 
                  oOctree::const_iterator children_beg, 
                  oOctree::const_iterator children_end) const override;

private:
  const oMapping& mapping;
  const oLevelSet& ls;
  const int dim;
  int levelMin, targetLevelMax;
};

/// Refine only cells with positive LS
class oRefinementCriteriaPositive : public oRefinementCriteria
{
public:
  oRefinementCriteriaPositive(const oMapping& m, const oLevelSet& ls_);
  bool operator()(oOctree::const_iterator cell, int level, 
                  const int * ijk, 
                  oOctree::const_iterator children_beg, 
                  oOctree::const_iterator children_end) const override;

private:
  const oMapping& mapping;
  const oLevelSet& ls;
  const int dim;
};


class oRefinementCriteriaOnBand : public oRefinementCriteria
{
public:
  oRefinementCriteriaOnBand(const oMapping& m, const oLevelSet& ls_, const double& d, int levelMin_ = 0);
  bool operator()(oOctree::const_iterator cell, int level, 
		  const int * ijk, 
		  oOctree::const_iterator children_beg, 
		  oOctree::const_iterator children_end) const override;
private:
  const oMapping& mapping;
  const oLevelSet& ls;
  const int dim;
  const double dist;
  const int levelMin;
};



class oRefinementCriteriaOnFunction : public oRefinementCriteria
{
public:
  oRefinementCriteriaOnFunction(const oMapping& m, const oLevelSet& ls_, std::function<bool (std::vector<double> &lsVals,
                                                                                             std::array<double, 3> &binf,
                                                                                             std::array<double, 3> &bsup)> func_, int levelMin_ = 0);
  bool operator()(oOctree::const_iterator cell, int level,
          const int * ijk,
          oOctree::const_iterator children_beg,
          oOctree::const_iterator children_end) const override;
private:
  const oMapping& mapping;
  const oLevelSet& ls;
  const int dim;
  const int levelMin;
  std::function<bool (std::vector<double> &lsVals,
                      std::array<double, 3> &binf,
                      std::array<double, 3> &bsup)> func;
};





class oRefinementCriteriaOnBandAndPositiveLevelset : public oRefinementCriteria
{
public:
  oRefinementCriteriaOnBandAndPositiveLevelset(const oMapping& m);
	void addLevelSet(const oLevelSet * ls, double d, int level, int sign=1);
  bool operator()(oOctree::const_iterator cell, int level, 
		  const int * ijk, 
		  oOctree::const_iterator children_beg, 
		  oOctree::const_iterator children_end) const override;

private:
  const oMapping& mapping;
  const int dim;
	std::vector<const oLevelSet* > ls_list;
	std::vector<double> d_list;
	std::vector<int> level_list;
	std::vector<int> sign_list;
};


//We want h< factor * curvature_radius
class oRefinementCriteriaOnCurvature : public oRefinementCriteria
{
public:
  oRefinementCriteriaOnCurvature(const oMapping& m, const oLevelSet& ls_, const double& fact);
  bool operator()(oOctree::const_iterator cell, int level,
                  const int * ijk,
                  oOctree::const_iterator children_beg,
                  oOctree::const_iterator children_end) const override;
private:
  const oMapping& mapping;
  const oLevelSet& ls;
  const int dim;
  const double factor;
};


} //end namespace 

#endif





