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
#ifndef H_SListBase
#define H_SListBase

struct SLinkBase{
  SLinkBase * next;
  SLinkBase() { next = nullptr; }
  SLinkBase(SLinkBase *p) { next = p; }
};

class SListBase{
protected:
  SLinkBase *last; //last->next is head
  int Size;
public:
  SListBase()
    : last(nullptr), Size(0) {} 
  SListBase(SLinkBase *a)
    : last(a->next=a), Size(1) {}
  
  int size() const { return Size; }
  void insert(SLinkBase *a);
  void append(SLinkBase *a);
  
  SLinkBase* get();
  
  void clear();
  SLinkBase* nth(int i) const;
  void reverse();
  
  ~SListBase();
  
  friend class SListBaseIter;
};

class SListBaseIter{
protected:
  SListBase *cs;
  SLinkBase *ce,*cl;
public:
  inline SListBaseIter(SListBase &s)
  { cs = &s; ce = cs->last; cl = nullptr;}
  inline void reset() { ce = cs->last; cl = nullptr; }
  void insert(SLinkBase *);
  inline void insertBefore(SLinkBase *);
  inline SLinkBase* operator() ();
  SLinkBase * next();
  void remove();
  SListBase const * list() const;
};

inline SLinkBase * SListBaseIter::operator() ()
{
  SLinkBase *ret = nullptr;
  if(ce)
    ret = (!(ce==cs->last && cl)) ? (cl=ce,ce=ce->next) : nullptr;
  return ret;
}

inline void SListBaseIter::insertBefore(SLinkBase *ne)
{
  if(cl){  // if there is a last element
    ne->next = ce;
    cl->next = ne;
    cs->Size++;
  } else // otherwise insert at start of list
    cs->insert(ne);
}

#endif
