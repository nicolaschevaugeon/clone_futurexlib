/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _VAL_MANAGER_H
#define _VAL_MANAGER_H

#include <cassert>
#include <algorithm>
#include <utility>
#include <functional>
#include <vector>
#include <map>
#include <string>
#include "eXlibris_hash_map.h"

#include <fstream>
#include <iostream>
#include "xDebug.h"


namespace xfem
{
template<typename KEY, typename VAL, typename HASH, typename EQUAL>
class xMapContainer {

public:
  typedef  eXlibris_types::hash_map<KEY, VAL*, HASH, EQUAL> container_t;
  typedef typename container_t::iterator       iterator;
  typedef typename container_t::const_iterator const_iterator;
  typedef typename container_t::value_type     value_type;

  container_t vals;
 
public:
 
  iterator begin()             { return vals.begin();}
  iterator   end()             { return vals.end();}
  const_iterator begin() const { return vals.begin();}
  const_iterator   end() const { return vals.end();}
  int             size() const { return vals.size();}

//old version to be removed
  //void Add(const KEY& d, VAL* v) {vals.insert(std::make_pair(d, v));}


  std::pair<iterator, bool> insert(const value_type& v) 
    {
       return vals.insert(v);
    }

  iterator find(const KEY& key) { return vals.find(key);}

  VAL* findSecond(const KEY& key) {
     iterator it = vals.find(key);
     if (it == vals.end()) return nullptr;
    return it->second;
  }

  const KEY* findFirst(const KEY& key) const
  {
     const_iterator it = vals.find(key);
     if (it == vals.end()) return nullptr;
    return &it->first;
  }
  void   erase (const KEY& key)  {
    iterator it = vals.find(key);
    if (it == vals.end()) return;
    delete it->second;
    vals.erase(it);
  }

  void clear() {
    iterator it = vals.begin(), itEnd = vals.end();
    for ( ; it != itEnd; ++it) delete it->second;
    vals.clear();
  }

  ~xMapContainer()   {
    iterator it = vals.begin(), itEnd = vals.end();
    for ( ; it != itEnd; ++it) delete it->second;
  }
};


//rule is pair Key, xValue* is in the map,
//   xValue* cannot be a null pointer
template <typename K, typename V, typename HASH, typename EQUAL>
class xValManager {

//typdefs

  //gestion des subsets
  //so far only one subset is coded
private:
  typedef typename std::map<std::string, std::vector<V*> *> subsets_t;
  typedef typename subsets_t::const_iterator const_iter_subsets;
  typedef typename subsets_t::iterator       iter_subsets;
  typedef typename std::vector<V*> subset_t;
  typedef xMapContainer<K, V, HASH, EQUAL> MapContainer_t;
  
public :
  typedef typename subset_t::const_iterator  vIter_const;
  typedef typename subset_t::iterator        vIter;
  typedef typename MapContainer_t::const_iterator map_const_iterator;
  typedef typename MapContainer_t::iterator map_iterator;
  typedef  V Value;



  ~xValManager()
    {
      MapContainer.clear();
      for (iter_subsets it = subsets.begin(); it != subsets.end(); ++it) 
	delete it->second;
    }

	V * getValPtr(const K& key) 
	{
		V* p = MapContainer.findSecond(key);
		if (!p) 
		{
			std::cerr << "no values associated to the key" << std::endl;
			std::cerr << key << std::endl;
			assert(1 == 0);
		}
		return p;
	}

  V * find(const K& key) 
    {
      V* p = MapContainer.findSecond(key);
      return p;
    }
  const K * findKeyAdress(const K& key) const 
  {
      const K* p = MapContainer.findFirst(key);
      return p;
  }

  template <class iterKey, class F>
  void insert(iterKey first, iterKey last, F& fct)
    {
      for (; first != last; ++first) insert(*first, fct);
    }

  template <class iterKey, class F>
  void insert(iterKey first, iterKey last, F& fct, std::vector<V*>& res)
    {
      for (; first != last; ++first) res.push_back(insert(*first, fct));
    }

  //if the key already exists, no insertion is performed
  //and the existing pointer is returned
  template <class F>
  V* insert(const K& key, F& fct)
    {
      map_const_iterator it = MapContainer.find(key);
      if (it != MapContainer.end()) return it->second;
      V* v = fct(key);
      if (v == nullptr) return nullptr;
//cout << "before val container insert" << endl; 
     MapContainer.insert(std::make_pair(key, v));
//cout << "after val container insert" << endl; 
//cout << "in insert return value is " << (p.first->second != 0) << endl;
//cout << "insertion was done  " << (p.second) << endl;
      return v;
    }

  void insert(const K& key, V* v)
    {
      MapContainer.insert(std::make_pair(key, v));
    }

  //careful !! the function does not remove the value in the subset
  //it should be coded
  void erase     (const K& key)        { MapContainer.erase(key); }

  void clear ()                    
  { 
    MapContainer.clear();  
    for (iter_subsets it = subsets.begin(); it != subsets.end(); ++it) 
      delete it->second;
    subsets.clear();
  }

  template <class iterKey, class T>
  void getValPtr(iterKey first, iterKey last, T& result) 
    {
      for (; first != last; ++first) result.push_back(getValPtr(*first));
    } 
  //io
  void PrintForDebug(const std::string& filename)
    {
      const bool debug = xdebug_flag;
      if (debug) std::cout << "taille MapContainer " << MapContainer.size() 
                           << std::endl;
      
      std::ofstream  out(filename.c_str());
      PrintForDebug(out);
    }

  template <class S, class Op> void getAll(std::vector<std::pair<K,S> > &buf)
    {     
      map_const_iterator it = MapContainer.begin(), itEnd = MapContainer.end();
      Op oper;
      for ( ; it != itEnd; ++it) 
      {
	typename Op::result_type r=oper(it->second);
	if (S dof=dynamic_cast<S>(r))
	{
	  buf.push_back(std::pair<K,S>(it->first,dof));
	}
      }
    }

  void PrintForDebug(std::ostream& out)
    {
      map_const_iterator it = MapContainer.begin(), itEnd = MapContainer.end();
      for ( ; it != itEnd; ++it) {
	out << it->first << std::endl;
	it->second->print(out) << std::endl;
      }
    }

   
   map_iterator begin()  { return MapContainer.begin(); }
   map_iterator  end()   { return MapContainer.end(); }
   map_const_iterator begin() const  { return MapContainer.begin(); }
   map_const_iterator  end() const   { return MapContainer.end(); }


// to manage subsets
   vIter begin(const std::string& sub)  { return getSubset(sub)->begin(); }
   vIter  end(const std::string& sub)   { return getSubset(sub)->end(); }
   vIter_const begin(const std::string& sub) const { return getSubset(sub)->begin(); }
   vIter_const  end(const std::string& sub) const   { return getSubset(sub)->end(); }
   int    size(const std::string& sub) const { return getSubset(sub)->size(); }
   void    add(V* v, const std::string& sub) {getSubset(sub)->push_back(v); }

/*   void    clear( const std::string& sub) {
     subset_t* s = getSubset(sub);
     for (typename std::vector<V*>::iterator it = s->begin(); it != s->end(); it++){
       (*it)->delState();
     }
     getSubset(sub)->clear();
   }*/
   void    clear_subset( const std::string& sub) {
     subset_t* s = getSubset(sub);
     delete s;
     subsets.erase(sub);
   }
   void    clear_all_subsets() {
     iter_subsets it = subsets.begin();
     iter_subsets ite = subsets.begin();
     for (; it != ite; ++it)
       { delete it->second; }
     subsets.clear();
   }

   void    erase( V* v, const std::string& sub) 
   {
        subset_t* s = getSubset(sub);
        vIter start = s->begin(); 
        vIter last = s->end(); 
        vIter f   = std::find (start,last,v); 
        if (f!=last)
        {
          s->erase(f);
        }
   }
  template <class F>
  V* update(const K& key, F& fct)
    {
      map_iterator it = MapContainer.find(key);
      V* v = fct(key);
      if (v == nullptr) return nullptr;
      if (it != MapContainer.end()) 
      {
          if (it->second != v)
          {
             delete it->second;
             it->second=v;
          }
      }
      else
      {
          MapContainer.insert(std::make_pair(key, v));
      }
      return v;
    }

  template <class F>
  xValManager& clone (const xValManager &other, const F& fct)
      {
          // clear every things
          clear(); 

          // loop to duplicate values of other into curent instance
          std::map<V *,V *> coresp;
          std::list<typename MapContainer_t::value_type> remains;
          for (auto v_pair : other.MapContainer )
          {
              V* new_val= fct (v_pair.second,coresp);
              if (new_val)
              {
                  coresp.insert(std::make_pair(v_pair.second,new_val));
                  MapContainer.insert(std::make_pair(v_pair.first, new_val));
              }
              else
                  remains.push_back(v_pair);
          }
          // finishing duplication for interdepandante values
          int old=remains.size();
          while (old)
          {
              auto it  = remains.begin();
              auto ite = remains.end();
              for (;it!=ite;)
              {
                  auto v_pair = *it;
                  V* new_val= fct (v_pair.second,coresp);
                  if (new_val)
                  {
                      coresp.insert(std::make_pair(v_pair.second,new_val));
                      MapContainer.insert(std::make_pair(v_pair.first, new_val));
                      remains.erase(it++);
                  }
                  else
                      ++it;
              }
              if (old==remains.size())
              {
                  std::cout << "Interdepandante value can not be duplicated " <<std::endl;
                  throw -24533;
              }
              else
                  old=remains.size();
          }
          // loop on subsets
          iter_subsets it = other.subsets.begin();
          iter_subsets ite = other.subsets.end();
          for (; it != ite; ++it)
          {
              subset_t* other_sub =  it->second;
              // create equivalent subset
              subset_t* sub = new subset_t;
              subsets.insert(std::make_pair(it->first, sub));

              // loop on value of other subset to fill new equivalent subset
              for (auto v : *other_sub)
              {
                  auto coresp_it = coresp.find(v);
                  assert(coresp_it!=coresp.end());
                  sub->push_back(coresp_it->second);
              }
          }
          return *this;
      }
private :  
  MapContainer_t MapContainer;
  mutable subsets_t subsets;
  subset_t* getSubset(const std::string& sub) const
    {
      const bool debug = xdebug_flag;
      const_iter_subsets it = subsets.find(sub);
      if (it == subsets.end()) {
	if (debug) std::cerr << "The entities set with name " 
                             << sub  << " does not exist\n";
	subset_t* ptr = new subset_t;
	subsets.insert(std::make_pair(sub, ptr));
	return ptr;
      }
      return it->second;
    }
};

} // end of namespace

#endif









