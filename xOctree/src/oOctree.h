/*
  octree is a subproject of  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

// -*- C++ -*-

#ifndef _OOCTREE_H__
#define _OOCTREE_H__

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>  //memcpy
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

// xoctree
#include "oMapping.h"

namespace xoctree
{
class oFilterCriteria;
class oRefinementCriteria;
class oObserver;
class oObserverRecorder;

/// Operation on ribbon to compare two states of the octree and return only the common active cells
class compareStatesOperation : public std::binary_function<int, int, int>
{
  public:
   int operator()(const int& f, const int& s) const { return (((f > 0) && (f != s)) * f); }  // New
};

/// Linear octree class
class oOctree
{
  public:
   typedef unsigned char cell_type;
   typedef cell_type* iterator;
   typedef const cell_type* const_iterator;

  public:
   oOctree(const oMapping& m_, int _LevelMax, int* per = nullptr, bool construct_ = true);
   virtual ~oOctree();

   const oMapping& getMapping() const { return mapping; }

   int getLevelMax() const { return LevelMax; }
   int getNbChildren() const { return NbChildren; }
   int getDim() const { return mapping.getDim(); }
   /// Get the periodicity vector
   const int* getPeriodicity() const { return Periodicity; }
   /// Number of cells for a given level
   inline int size(int level) const { return nsize[level]; }
   inline iterator begin(int level)
   {
      assert(level <= LevelMax);
      return cells + ibegin[level];
   }
   inline iterator end(int level)
   {
      assert(level <= LevelMax);
      return cells + iend[level];
   }
   inline const_iterator begin(int level) const
   {
      assert(level <= LevelMax);
      return cells + ibegin[level];
   }
   inline const_iterator end(int level) const
   {
      assert(level <= LevelMax);
      return cells + iend[level];
   }
   /// Get the total size of the octree (active and non active cells)
   inline int size() const { return CumulativeSize; }
   inline iterator begin() { return cells; }
   inline iterator end() { return cells + iend[LevelMax]; }
   inline const_iterator begin() const { return cells; }
   inline const_iterator end() const { return cells + iend[LevelMax]; }

   /// Print the activity of the octree (state of the cells from level 0 to level max)
   void printActivity() const;

   /// Compress the octree storage
   template <class ITER>
   static void compress(ITER it, ITER ite, std::set<unsigned int>& compressed)
   {
      const bool debug = false;
      ITER itb = it;
      compressed.clear();
      int count = 0;
      for (; it != ite; ++it)
      {
         if (*it == 1) compressed.insert(compressed.end(), count);
         count++;
      }
      if (debug)
      {
         std::cout << " activity " << std::endl;
         std::copy(itb, ite, std::ostream_iterator<unsigned int>(std::cout, " "));
         std::cout << std::endl << " compressed activity " << std::endl;
         std::copy(compressed.begin(), compressed.end(), std::ostream_iterator<int>(std::cout, " "));
         std::cout << std::endl;
      }
      std::cout << " compression efficiency " << (double)count / (double)compressed.size() << std::endl;
   }

   void compress_descendants(const_iterator cell, int level, std::set<unsigned int>& compressed) const;

   /// \brief get the cartesian coords (ijk) ON THE FINEST LEVEL of a given cell with given level.
   /// \param[in] cell : cell pointer
   /// \param[in] level : level of the cell
   /// \param[out] ijk_fine : int[3] on finest level
   inline void octree2finest_cartesian(const_iterator cell, int level, int* ijk_fine) const
   {
      while (level != LevelMax) cell = getChildren(cell, level++);
      octree2cartesian(cell, LevelMax, ijk_fine);
   }

   /// Filter the cells (keep only active cells that fullfill the filter criterion)
   void filter_cells(oFilterCriteria& filter_criteria);
   /// optimize the octree shape with respect to a given criterion (FROM SCRATCH)
   void optimize(oRefinementCriteria& test, bool twoForOneRule = true);
   /// refine the octree with respect to a given criterion
   void refine(oRefinementCriteria& refinement_criteria);
   /// refine the octree with respect to a given criterion (only for cells between level ldeb and lfin)
   void refine(oRefinementCriteria& refinement_criteria, int ldeb, int lfin);

   /// Get childrens of a cell.
   /// \param[in] cell : cell
   /// \param[in] level : level of the sons (of the cell ??)
   inline cell_type* getChildren(const_iterator cell, int level) const
   {
      const bool debug = false;
      if (debug) std::cout << " get children for level " << level << std::endl;
      if (debug) std::cout << " for parent location " << cell - begin(level) << std::endl;
      if (debug)
         std::cout << " children location in the global activity vector are "
                   << begin(level + 1) + NbChildren * (cell - begin(level)) - cells << std::endl;
      return (cell_type*)begin(level + 1) + NbChildren * (cell - begin(level));
   }

   /// Get Number of Active children of a cell. (all active childrens from level +1  to level_max)
   /// \param[in] cell : cell
   /// \param[in] level : level of the cell
   /// \param[in] level_max : max lineage level
   /// \param[in] generation: for recursivity
   inline int getNumberOfActiveChidren(const_iterator cell, int level, int level_max, int generation = 0)
   {
      int sons = 0;
      // If cell is already active at first generation, then nothing to do
      if (generation == 0 && *cell)
      {
         return 0;
      }

      cell_type* childrens = getChildren(cell, level);
      for (int ic = 0; ic < NbChildren; ++ic)
      {
         cell_type* child = childrens;
         if (*child)
            sons += 1;
         else if (level <= LevelMax)
            sons += getNumberOfActiveChidren(childrens, level + 1, level_max, generation + 1);
         ++childrens;
      }
      return sons;
   }

   /// Get Active Lineage of a cell. (all active childrens from level +1  to level_max)
   /// \param[in] cell : cell
   /// \param[in] level : level of the cell
   /// \param[in] level_max : max lineage level
   /// \param[out] lineage: list of pair<cell*,int> of active sons
   /// \param[in] generation: for recursivity
   inline void getActiveLineage(const_iterator cell, int level, int level_max, std::list<std::pair<cell_type*, int>>& lineage,
                                int generation = 0)
   {
      // If cell is already active at first generation, then nothing to do
      if (generation == 0 && *cell)
      {
         //          lineage.push_back(std::make_pair(cell,level));
         return;
      }

      cell_type* childrens = getChildren(cell, level);
      for (int ic = 0; ic < NbChildren; ++ic)
      {
         cell_type* child = childrens;
         if (*child)
            lineage.push_back(std::make_pair(child, level + 1));
         else if (level <= LevelMax)
            getActiveLineage(childrens, level + 1, level_max, lineage, generation + 1);
         ++childrens;
      }
      return;
   }

   /// Get ancestor for a cell
   /// \param[in] up : level difference between cell and ancestor
   cell_type* getAncestor(int up, const cell_type* cell, int level);

   //  void getNeighbors  (const int* ijk, int level, int* neighbors,
   //		      int& nb_neighbors );

   /// Loop on cells of a given level
   /// \param[out] ijk : ijk of the next cell
   /// \param[in] level : level of interest
   /// \return false after last cell
   bool next_ijk(int* ijk, int level) const;

   /// Is the cell on the boundary ?
   /// \param[in] b : number of the bnd (1 to 4 for 2D, 1 to 6 for 3D)
   /// \param[in] ijk : ijk of the cell
   /// \param[in] level : level of the cell
   bool is_on_boundary(int b, const int* ijk, int level) const;

   /// Activate a cell in the database
   void activate(cell_type* cell, int level, const int* ijk);
   /// Deactivate a cell in the database
   void deactivate(cell_type* cell, int level, const int* ijk);
   /// Deactivate finest level
   void deactivateFinestLevel();
   void activateFinestLevel();

   /// Translates octree coords of the cell (pointer + level) to cartesian coords (ijk + level)
   virtual void octree2cartesian(const_iterator cell_octree, int level, int* ijk) const;
   /// Translates cartesian coords of the cell (ijk + level) to octree coords (pointer + level)
   virtual cell_type* cartesian2octree(const int* ijk, int level) const;

   /// How many cells are active ?
   int getNbActiveCells() const;

   /// Deprecated ????
   void refineOctree(oRefinementCriteria& refinement_criteria);

   /// Euhhhh... activate a cell, but not her sons ??? (je ne me rappelle plus !!!)
   void activateComp(oOctree::cell_type* cell, int level, const int* ijk);

   /// Compare two states of the octree based on the record of the modifications, and activate ONLY what was modified.
   void compareOctree(oObserverRecorder& obs);

   // Need to be reimplemented !

   void saveState()
   {
      //    cellsSaved = new cell_type [ CumulativeSize ];
      //    memcpy(cellsSaved,cells,CumulativeSize);
      cellsSaved.resize(CumulativeSize);
      std::copy(cells, cells + CumulativeSize, cellsSaved.begin());
   }

   void restoreState()
   {
      //    memcpy(cells,cellsSaved,CumulativeSize);
      //    delete cellsSaved;
      std::copy(cellsSaved.begin(), cellsSaved.end(), cells);
      cellsSaved.resize(0);
   }

#if 0
    void saveState(std::string name){
      savedStates[name]= (new cell_type [ CumulativeSize ]);
      memcpy(savedStates[name],cells,CumulativeSize);
    }
  
    void restoreState(std::string name){
      memcpy(cells,savedStates[name],CumulativeSize);
      delete savedStates[name];/*Erase memory*/
      //Erase map entry :
      std::map<std::string, cell_type *>::iterator iter = savedStates.find(name);
      savedStates.erase(iter);
    }
  
    void eraseState(std::string name){
      delete savedStates[name];/*Erase memory*/
      //Erase map entry :
      std::map<std::string, cell_type *>::iterator iter = savedStates.find(name);
      savedStates.erase(iter);
    }
    
    void compareStates(std::string name1, std::string name2="current"){
      // transform(begin(), end(), savedStates[name1], begin(), std::minus<int>());
      // transform(begin(), end(), cellsSaved, begin(), std::minus<int>());
      // transform(begin(), end(), cellsSaved, begin(), compareStates_c());
    
      if(name2=="current"){
	transform(begin(), end(), savedStates[name1], begin(), compareStatesOperation());
      }else{
	transform(savedStates[name2], savedStates[name2]+iend[LevelMax], 
		  savedStates[name1], begin(), compareStatesOperation());
      }
    }
#endif

   // void getNodeIdsOnElt (int* ijk_ori,  int level, int *  node_ids) const;

   /// Print the octree vector
   void viewState();

   /// Attach an observer
   void attach(oObserver& obs_);
   /// Detach an observer
   void detach(oObserver& obs_);
   /// notify refinement to observers
   void notify_refinement(const cell_type* cell, int level, const int* ijk);
   /// notify derefinement to observers
   void notify_derefinement(const cell_type* cell, int level, const int* ijk);

   void lookForNodeOnNeighbors2D(const cell_type* cell, int level, const int* ijk, int* exixt_on_edge) const;
   void lookForNodeOnNeighbors3D(const cell_type* cell, int level, const int* ijk, int* exixt_on_edge, int* exixt_on_face) const;
   void lookForNodeHangingOnNeighbors2D(const cell_type* cell, int level, const int* ijk, int* hanging_on_edge) const;
   void lookForNodeHangingOnNeighbors3D(const cell_type* cell, int level, const int* ijk, int* hanging_on_edge,
                                        int* hanging_on_face) const;
   void refineForOneLevelDiffConstraint();

   /// Put octree back to its original state (modifications recorded by the observer)
   void modifyOctree(oObserverRecorder& obs);
   const oTopo& getTopo() const { return topo; }

   /// Reinit octree shape (initial levels along x, y and z)
   /// Warning : the octree structure is completely resetted (both state AND pointers !!), use reset if this is unwanted !
   //! warning : with oOctreeAnisotropic derived class, using here some new motif (different from object creation) looks
   //! dangerous as this methode doesn't modify the private member call motif.
   void reInitOctree(const int* motif = nullptr) { construct(motif); }
   void reset();

   void locatePointOnFinestLevel(double x, double y, double z, int* ijk);
   void locatePointOnActiveCells(double x, double y, double z, int* ijk, int& level);

  protected:
   /// Refine the octree without taking care of the 2:1 rule (between level ldeb and lfin)
   void refine_incompatible(oRefinementCriteria& test, int ldeb, int lfin);
   /// Coarsen the octree without taking care of the 2:1 rule
   void coarsen_incompatible(oRefinementCriteria& refinement_criteria);
   /// Construct the database
   void construct(const int* motif = nullptr);

   const oMapping& mapping;
   const int LevelMax;
   const int Dim;
   int Periodicity[3];
   const int NbChildren;
   int CumulativeSize;
   cell_type* cells;
   //   int nsize[topo.MAX_DEPTH];
   //   int ibegin[topo.MAX_DEPTH];
   //   int iend[topo.MAX_DEPTH];
   int* nsize;
   int* ibegin;
   int* iend;
   int* ijk2_almost_octree;
   //  cell_type *cellsSaved;
   std::vector<cell_type> cellsSaved;
   std::map<std::string, cell_type*> savedStates;
   std::minus<int> moins;
   std::list<oObserver*> obs;
   const oTopo& topo;
};

/// Criteria base class to filter the cells of the octree
class oFilterCriteria
{
  public:
   virtual ~oFilterCriteria() = default;
   virtual bool operator()(oOctree::const_iterator cell, int level, const int* ijk, oOctree::const_iterator children_beg,
                           oOctree::const_iterator children_end) const = 0;

  private:
};

/// Anisotropic octree
class oOctreeAnisotropic : public oOctree
{
  public:
   oOctreeAnisotropic(const oMapping& m_, int _LevelMax, int* per = nullptr);
   ~oOctreeAnisotropic() override;

   void octree2cartesian(const_iterator cell_octree, int level, int* ijk) const override;
   cell_type* cartesian2octree(const int* ijk, int level) const override;

   //  void locatePointOnFinestLevel(double x, double y, double z, int *ijk);
   //  void locatePointOnActiveCells(double x, double y, double z, int *ijk, int &level);

  private:
   const int* motif;
   bool anisotropic;

   // Constantes utiles :
   const int layer1D_0;         // nb de cellules selon x au niveau 0
   const int layer2D_0;         // nb de cellules selon (x,y) au niveau 0
   const int* generation_size;  // pointer to vector of value depending only on level. Dim is treated at construction
};

/// Criteria base class to refine the cells of the octree
class oRefinementCriteria
{
  public:
   virtual ~oRefinementCriteria() = default;
   virtual bool operator()(oOctree::const_iterator cell, int level, const int* ijk, oOctree::const_iterator children_beg,
                           oOctree::const_iterator children_end) const = 0;

  private:
};

class oRefinementCriteriaList : public oRefinementCriteria
{
  public:
   void insert(const oRefinementCriteria& r) { criterialist.push_back(&r); };
   void clear() { criterialist.clear(); };
   bool operator()(oOctree::const_iterator cell, int level, const int* ijk, oOctree::const_iterator children_beg,
                   oOctree::const_iterator children_end) const override
   {
      bool returnval = false;
      std::list<const oRefinementCriteria*>::const_iterator it = criterialist.begin();
      while ((it != criterialist.end()) && !returnval)
      {
         returnval += (*it)->operator()(cell, level, ijk, children_beg, children_end);
         ++it;
      }
      return returnval;
   };

  private:
   std::list<const oRefinementCriteria*> criterialist;
};

class oRefinementCriteriaComplement : public oRefinementCriteria
{
  public:
   oRefinementCriteriaComplement(oRefinementCriteria& criterion_) : criterion(criterion_){};
   bool operator()(oOctree::const_iterator cell, int level, const int* ijk, oOctree::const_iterator children_beg,
                   oOctree::const_iterator children_end) const override
   {
      return (!criterion(cell, level, ijk, children_beg, children_end));
   };

  private:
   oRefinementCriteria& criterion;
};

class oRefinementCriteriaOnLevel : public oRefinementCriteria
{
  public:
   oRefinementCriteriaOnLevel(int targetLevel_) : targetLevel(targetLevel_){};
   bool operator()(oOctree::const_iterator cell, int level, const int* ijk, oOctree::const_iterator children_beg,
                   oOctree::const_iterator children_end) const override
   {
      return (level < targetLevel);
   };

  private:
   int targetLevel;
};

class oRefinementCriteriaTrue : public oRefinementCriteria
{
  public:
   oRefinementCriteriaTrue(){};
   bool operator()(oOctree::const_iterator cell, int level, const int* ijk, oOctree::const_iterator children_beg,
                   oOctree::const_iterator children_end) const override
   {
      return (true);
   };

  private:
};

/// Observer base class : used to monitor the modifications of the octree
class oObserver
{
  public:
   virtual ~oObserver() = default;
   virtual void refine(const oOctree& octree, const oOctree::cell_type* cell, int level, const int* ijk) = 0;
   virtual void derefine(const oOctree& octree, const oOctree::cell_type* cell, int level, const int* ijk) = 0;

  private:
};

class oObserverRecorder : public oObserver
{
  public:
   oObserverRecorder() {}
   void refine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk) override;
   void derefine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk) override;
   void clear() { record.clear(); }
   int size() { return record.size(); }

   //   bool findSons(oOctree::cell_type *cell, int level);
   //   void fillSons();
   void getSons(const oOctree& octree, std::set<std::vector<int>>& sons);
   void getGrandSons(const oOctree& octree, std::set<std::vector<int>>& gsons);

   typedef std::vector<std::pair<std::pair<int, std::vector<int>>, bool>> record_type;
   record_type::iterator begin() { return record.begin(); }
   record_type::iterator end() { return record.end(); }

  private:
   record_type record;
   //   record_type recordSons;
};

//   template<class ITER>
//   void oOctree::modifyOctree(ITER it, ITER end){
//   const bool debug=true;
//     for(;it!=end;++it) {
//       bool state= it->second;
//       int level= (it->first).first;
//       std::vector<int> ijk=(it->first).second;
//       if(debug){
// 	cout<<"state="<<state<<" level="<<level<<" ijk="<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<std::endl;
//       }
//       cell_type* cell=cartesian2octree(&ijk[0],level);
//       if(state) activate(cell,level,&ijk[0]);
//       else deactivate(cell,level,&ijk[0]);
//       }
//   }

}  // namespace xoctree

#endif
