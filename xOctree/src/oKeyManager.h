/*
  octree is a subproject of  xfem : C++ Finite Element Library
  developed under the GNU Lesser General Public License
  See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

// -*- C++ -*-

#ifndef _OKEYMANAGER_H__
#define _OKEYMANAGER_H__

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>

#include "oOctree.h"

namespace xoctree
{
class oKey
{
  public:
   oKey(const int* ijk_) : is_regular(true), is_fixed(false) { std::copy(ijk_, ijk_ + 3, ijk); }
   bool operator<(const oKey& other) const
   {
      if (ijk[0] < other.ijk[0]) return true;
      if (ijk[0] > other.ijk[0]) return false;
      if (ijk[1] < other.ijk[1]) return true;
      if (ijk[1] > other.ijk[1]) return false;
      if (ijk[2] < other.ijk[2]) return true;
      return false;
   }
   const int* getIJK() const { return ijk; }

   typedef std::vector<oKey*>::const_iterator const_iterator;
   typedef std::vector<oKey*>::iterator iterator;

   const_iterator beginTies() const
   {
      assert(!is_regular && (keys.size() > 1));
      return keys.begin();
   }
   iterator beginTies()
   {
      assert(!is_regular && (keys.size() > 1));
      return keys.begin();
   }
   const_iterator endTies() const
   {
      assert(!is_regular && (keys.size() > 1));
      return keys.end();
   }
   iterator endTies()
   {
      assert(!is_regular && (keys.size() > 1));
      return keys.end();
   }
   int sizeTies() const
   {
      assert(!is_regular && (keys.size() > 1));
      return keys.size();
   }
   const std::vector<oKey*>& getTies() const
   {
      assert(!is_regular && (keys.size() > 1));
      return keys;
   }
   const oKey* getMirror() const
   {
      assert(!is_regular && (keys.size() == 1));
      return keys[0];
   }
   oKey* getMirror()
   {
      assert(!is_regular && (keys.size() == 1));
      return keys[0];
   }

   const_iterator beginStencil() const
   {
      assert(is_regular);
      return keys.begin();
   }
   iterator beginStencil()
   {
      assert(is_regular);
      return keys.begin();
   }
   const_iterator endStencil() const
   {
      assert(is_regular);
      return keys.end();
   }
   iterator endStencil()
   {
      assert(is_regular);
      return keys.end();
   }
   int sizeStencil() const
   {
      assert(is_regular);
      return keys.size();
   }
   const std::vector<oKey*>& getStencil() const
   {
      assert(is_regular);
      return keys;
   }

   friend std::ostream& operator<<(std::ostream& ofs, const oKey& k)
   {
      ofs << "ijk : ";
      std::copy(k.ijk, k.ijk + 3, std::ostream_iterator<int>(ofs, " "));
      ofs << " id is " << k.getId() << std::endl;
      if (k.is_regular)
         ofs << " the node is regular" << std::endl;
      else
      {
         ofs << " the node is hanging or periodic" << std::endl;
         ofs << " the ids of the related keys are \n  ";
         for (const_iterator it = k.keys.begin(); it != k.keys.end(); ++it) ofs << (*it)->getId() << " ";
         ofs << std::endl;
      }
      ofs << " Value is : " << std::endl;
      unsigned int j = k.vals.size();
      for (unsigned int i = 0; i < j; ++i)
      {
         oKey::setField(i);
         ofs << "     Field " << i << " : " << k.getVal() << std::endl;
      }
      return ofs;
   }
   void setHanging(oKey* k1, oKey* k2)
   {
      is_regular = false;
      keys.clear();
      keys.push_back(k1);
      keys.push_back(k2);
   }
   void setHanging(oKey* k1, oKey* k2, oKey* k3, oKey* k4)
   {
      is_regular = false;
      keys.clear();
      keys.push_back(k1);
      keys.push_back(k2);
      keys.push_back(k3);
      keys.push_back(k4);
   }
   void setRegular()
   {
      if (!is_regular) setValsFromTies();
      is_regular = true;
   }

   void setStencil(int i, oKey* k)
   {
      is_regular = true;
      keys[i] = k;
   }
   void resizeStencil(int n) { keys.resize(n, nullptr); }

   bool isRegular() const { return is_regular; }
   bool isPeriodic() const { return (!is_regular && (keys.size() == 1)); }
   bool isHanging() const { return (!is_regular && (keys.size() != 1)); }
   static double step(int i) { return step_[i]; }
   static void setCoeff(const int* c) { std::copy(c, c + 3, coeff_); }
   static void setStep(const double* c) { std::copy(c, c + 3, step_); }

   int getId() const { return std::inner_product(ijk, ijk + 3, coeff_, 0); }

   void setDof(int n) { NumDof = n; }
   bool isDof() const { return (is_regular && !is_fixed); }
   int getNumDof() const { return NumDof; }

   double getVal() const
   {
      const bool debug = false;
      if (is_regular)
      {
         if (debug) std::cout << " i_field is " << i_field << std::endl;
         return vals[oKey::i_field];
      }
      // one line for below ?
      double v = 0.;
      std::vector<oKey*>::const_iterator it = keys.begin();
      std::vector<oKey*>::const_iterator ite = keys.end();
      for (; it != ite; ++it) v += (*it)->getVal();
      return v / (double)keys.size();
   }
   std::vector<double>::const_iterator beginVal() const { return vals.begin(); }
   std::vector<double>::const_iterator endVal() const { return vals.end(); }
   std::vector<double>::iterator beginVal() { return vals.begin(); }
   std::vector<double>::iterator endVal() { return vals.end(); }

   void setVal(const double& v) { vals[oKey::i_field] = v; }
   void setVal(const std::vector<double>& v) { vals = v; }
   void pushVal(const double& v) { vals.push_back(v); }

   friend struct oEqualKey;
   friend double distance(const oKey* k1, const oKey* k2, int i);
   friend class oKeyManager;
   static void setField(int i) { oKey::i_field = i; }

  private:
   oKey() {}
   int ijk[3];
   bool is_regular;
   // if is_regular is true the stencil is stored in the "keys" member
   //   (pointer zero if no point in some direction)
   // if is_regular is false, the linear combination is stored in the "keys" member
   std::vector<oKey*> keys;  // store stencil or ties
   static double step_[3];
   static int coeff_[3];

   bool is_fixed;
   int NumDof;

   static int i_field;
   std::vector<double> vals;
   void setValsZero(int n) { vals.resize(n, 0.); }
   void getVals(std::vector<double>& vec) const
   {
      for (unsigned int i = 0; i < vals.size(); ++i)
      {
         oKey::setField(i);
         vec[i] = getVal();
      }
   }
   void setValsFromTies()
   {
      for (unsigned int i = 0; i < vals.size(); ++i)
      {
         oKey::setField(i);
         vals[i] = getVal();
      }
   }
};

struct oEqualKey
{
   bool operator()(const oKey* k1, const oKey* k2) const
   {
      return (k1->ijk[0] == k2->ijk[0] && k1->ijk[1] == k2->ijk[1] && k1->ijk[2] == k2->ijk[2]);
   }
};

struct oHashKey
{
   size_t operator()(const oKey* k) const { return k->getId(); }
};

double distance(const oKey* k1, const oKey* k2, int i);
void distance(const oKey* k, const std::vector<oKey*>& keys, std::vector<double>& hs);

struct oDeleteObject
{
   template <typename T>
   void operator()(const T* ptr) const
   {
      delete ptr;
   }
};

class oKeyManager
{
  public:
   typedef std::unordered_set<oKey*, oHashKey, oEqualKey> container;
   typedef container::const_iterator const_iterator;
   typedef container::iterator iterator;
   typedef container::size_type size_type;
   typedef container::value_type value_type;

  public:
   oKeyManager(const oOctree& octree_)
       : octree(octree_), dim(octree_.getDim()), level_max(octree_.getLevelMax()), nb_field(0), topo(octree.getTopo())
   {
      int c[3] = {1, (topo.index_max_for_levels[0][level_max] + 2),
                  (topo.index_max_for_levels[0][level_max] + 2) * (topo.index_max_for_levels[1][level_max] + 2)};
      oKey::setCoeff(c);
      const oMapping& mapping = octree.getMapping();
      oKey::setStep(mapping.getStep(level_max));
      createKeys();
   }
   ~oKeyManager()
   {
      std::for_each(keys.begin(), keys.end(), oDeleteObject());
      keys.clear();
   }
   void clear()
   {
      std::for_each(keys.begin(), keys.end(), oDeleteObject());
      keys.clear();
   }
   iterator begin() { return keys.begin(); }
   iterator end() { return keys.end(); }
   const_iterator begin() const { return keys.begin(); }
   const_iterator end() const { return keys.end(); }
   size_type size() const { return keys.size(); }

   std::pair<oKey*, bool> insert(const int* ijk)
   {
      oKey* key = new oKey(ijk);
      std::pair<iterator, bool> check = keys.insert(key);
      if (!check.second) delete key;
      key = *check.first;
      key->setValsZero(nb_field);
      return std::make_pair(key, check.second);
   }

   size_t erase(const int* ijk)
   {
      std::copy(ijk, ijk + 3, dummy.ijk);
      iterator it = keys.find(&dummy);
      delete (*it);
      keys.erase(it);
      return 1;
   }

   void eraseAll()
   {
      //    iterator it = keys.begin();
      //    iterator ite = keys.end();
      //    for (;it != ite; ++it) keys.erase(*it);
      keys.erase(keys.begin(), keys.end());  // Yurk yurk yurk !!!
      throw;                                 // needs proper deletion of pointers *oKey
   }

   const oKey* find(const int* ijk) const
   {
      std::copy(ijk, ijk + 3, dummy.ijk);
      oKeyManager::const_iterator it = keys.find(&dummy);
      if (it == keys.end()) return nullptr;
      return *it;
   }

   oKey* find(const int* ijk)
   {
      std::copy(ijk, ijk + 3, dummy.ijk);
      oKeyManager::iterator it = keys.find(&dummy);
      if (it == keys.end()) return nullptr;
      return *it;
   }

   void getKeysOnElt(const int* ijk_ori, int level, std::vector<const oKey*>& elt_keys) const;
   void getKeysOnElt(const int* ijk_ori, int level, std::vector<oKey*>& elt_keys);
   void getNodeIdsOnElt(int* ijk_ori, int level, int* node_ids) const;
   int getNodeId(int* ijk_ori, int level) const;
   void createKeys();
   void createStencils();

   typedef std::function<double(const int*)> f_ijk_t;
   int registerField(f_ijk_t f)
   {
      iterator it = begin();
      iterator ite = end();
      for (; it != ite; ++it) (*it)->pushVal(f((*it)->getIJK()));
      // the line below does not work because pushVal takes a reference
      // std::for_each(begin(), end(), bind2nd(std::mem_fun(&oKey::pushVal), 0.));
      nb_field++;
      return nb_field - 1;
   }
   int getNbField() const { return nb_field; }

   void refine(oOctree::const_iterator cell, int level, const int* ijk);
   void derefine(oOctree::const_iterator cell, int level, const int* ijk);
   bool isOnBoundary(oKey* key) const;
   bool isOnBBoxEdges(oKey* key) const;
   void printForDebug(const std::string& filename);
   void printForDebug(std::ostream& out);

   float max_load_factor() { return keys.max_load_factor(); }
   void max_load_factor(float lf) { return keys.max_load_factor(lf); }

  private:
   const oOctree& octree;
   container keys;
   mutable oKey dummy;
   const int dim;
   const int level_max;
   int nb_field;
   const oTopo& topo;

   // int exist_on_neighbor_edge[12];
   // int exist_on_neighbor_face[6];
};

class oObserverKeyManager : public oObserver
{
  public:
   oObserverKeyManager(oKeyManager& k_);
   void refine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk);
   void derefine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk);

  private:
   oKeyManager& key_manager;
};

}  // namespace xoctree

#endif
