/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#ifndef H_SSList
#define H_SSList

#include "SListBase.h"

template<class T>
struct SSLink : public SLinkBase{
  T info;
  SSLink(T a) : info(a) {}
};

template<class T> class SSListCIter;
template<class T> class SSListIter;

template<class T>
class SSList : private SListBase{
friend class SSListCIter<T>;
friend class SSListIter<T>;
public:
  SSList() : SListBase() {}
  SSList(SSLink<T> *a) : SListBase(a) {}

  SSList(const SSList<T> &c);

  int size() const
    { return SListBase::size(); }
  void insert(const SSList<T> &al);
  void insert(const T& a)
    { SListBase::insert( new SSLink<T>(a)); }
  void append(const T& a)
    { SListBase::append( new SSLink<T>(a)); }
  void append(const SSList<T> &al);
  void remove(const T &a);
  int inList(const T& a) const;
  void appendUnique(const T& a);
  void appendUnique(const SSList<T> &al);
  T& nth(int i) const
    { return ((SSLink<T>*)(SListBase::nth(i)))->info;}
  void clear()
    { SListBase::clear(); }
  void reverse()
    { SListBase::reverse(); }
  
  void operator=(const SSList<T> &l);
};

template<class T>
class SSListCIter : private SListBaseIter {
public:
  SSListCIter(const SSList<T>& s) : SListBaseIter( (SSList<T>&)s) {}

  inline int operator() (T &el)
  { 
    SSLink<T> *lnk = (SSLink<T>*) SListBaseIter::operator() ();
    if(lnk)
      el = lnk->info;
    return (lnk!=nullptr);
  }
  int next(T &el);
  inline void reset()
    { SListBaseIter::reset(); }
  SSList<T> const * list() const;
};


template<class T>
class SSListIter : private SListBaseIter {
public:
  SSListIter(SSList<T>& s) : SListBaseIter(s) {}

  inline T* operator() ()
    { SSLink<T> *lnk = (SSLink<T>*) SListBaseIter::operator() ();
      return lnk ? &lnk->info : nullptr; }
  inline int operator() (T &el)
    { SSLink<T> *lnk = (SSLink<T>*) SListBaseIter::operator() ();
      if(lnk){
	el = lnk->info;
      }
      return (lnk!=nullptr);
    }
  //return lnk ? el == (el = lnk->info) : 0; }
  inline void reset()
    { SListBaseIter::reset(); }
  inline void insert(T ne)
    { SListBaseIter::insert( new SSLink<T>(ne) ); }
  inline void insertBefore(T ne)
    { SListBaseIter::insertBefore( new SSLink<T>(ne) ); }
  inline void remove()
    { SListBaseIter::remove(); }
  void replaceCurrent( T ne );
};

template<class T>
SSList<T>::SSList(const SSList<T> &c) : SListBase()
{
  SSLink<T> *ce = (SSLink<T> *)(c.last);
  
  while(ce){
    ce = ce ? (SSLink<T> *)(ce->next) : nullptr;
    append( ce->info );
    if(ce == (SSLink<T> *)(c.last))
      ce = nullptr;
  }
}


template<class T>
void SSList<T>::operator=(const SSList<T> &c)
{
  clear();
  SSLink<T> *ce = (SSLink<T> *)(c.last);
  
  while(ce){
    ce = ce ? (SSLink<T> *)(ce->next) : nullptr;
    append( ce->info );
    if(ce == (SSLink<T> *)(c.last))
      ce = nullptr;
  }
}

template<class T>
int SSList<T>::inList(const T& a) const
{
  SSLink<T> *ce = (SSLink<T> *)last;
  
  while(ce){
    ce = ce ? (SSLink<T> *)(ce->next) : nullptr;
    if (ce->info == a)
      return 1;
    if(ce == (SSLink<T> *)last)
      ce = nullptr;
  }
  return 0;
}

template<class T>
void SSList<T>::remove(const T& remItem)
{
  SSListIter<T> iter(*this);
  T a;
  while(iter(a)){
    if(a == remItem){
      iter.remove();
      break;
    }
  }
}

template<class T>
void SSList<T>::insert(const SSList<T> &al)
{
  SSListCIter<T> aIter(al);
  SSListIter<T> thisIter(*this);
  T item;
  while(aIter(item)){
    thisIter.insert(item);
    thisIter();
  }
}

template<class T>
void SSList<T>::append(const SSList<T> &al)
{
  SSLink<T> *ce = (SSLink<T> *)(al.last);
  
  while(ce){
    ce = ce ? (SSLink<T> *)(ce->next) : 0;
    append(ce->info);
    if(ce == (SSLink<T> *)(al.last))
      ce = 0;
  }
}

template<class T>
void SSList<T>::appendUnique(const T& a)
{
  if( ! inList(a) )
    append(a);
}

template<class T>
void SSList<T>::appendUnique(const SSList<T> &al)
{
  SSLink<T> *ce = (SSLink<T> *)(al.last);
  
  while(ce){
    ce = ce ? (SSLink<T> *)(ce->next) : nullptr;
    appendUnique(ce->info);
    if(ce == (SSLink<T> *)(al.last))
      ce = nullptr;
  }
}

template<class T>
int SSListCIter<T>::next(T &el)
{ 
  SSLink<T> *lnk = (SSLink<T>*)SListBaseIter::next();
  if(lnk){
    el = lnk->info;
    return 1;
  } else
    return 0;
}

 template<class T>
 void SSListIter<T>::replaceCurrent( T ne )
 { ((SSLink<T>*)ce)->info = ne; }

template<class T>
SSList<T> const * SSListCIter<T>::list() const
{ return (SSList<T> *)SListBaseIter::list(); }


#endif
