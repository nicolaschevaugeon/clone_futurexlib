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
#include "SListBase.h"

SListBase::~SListBase()
{
  SLinkBase *cur = nullptr;
  while((cur = SListBase::get()))
    delete cur;
}


void SListBase::insert(SLinkBase *a)
{
  if(last)
    a->next = last->next;
  else
    last = a;
  last->next = a;
  Size ++;
}

void SListBase::append(SLinkBase *a)
{
  if(last){
    a->next = last->next;
    last = last->next = a;
  } else
    last = a->next = a;
  Size ++;
}

SLinkBase * SListBase::get()
{
  if(last == nullptr)
    return nullptr;
  SLinkBase *f = last->next;
  if(f == last)
    last = nullptr;
  else
    last->next = f->next;
  Size--;
  return f;
}

SLinkBase * SListBase::nth(int i) const
{
  SLinkBase *ce = last->next;
  int n = 0;
  while(ce){
    if(n==i)
      return ce;
    ce = ce ? ce->next : nullptr;
    if(ce == last->next)
      ce = nullptr;
    n++;
  }
  return nullptr;
}

void SListBase::reverse()
{
  SListBase temp;
  SLinkBase *ce;
  while((ce = get())){
    temp.insert(ce);
  }
  last = temp.last;
  Size = temp.Size;
  temp.last = nullptr;
  temp.Size = 0;
}


void SListBaseIter::remove()
{
  if(cl){
    if(cs->Size==1){ // deleting only item in list
      delete ce;
      cs->last = ce = cl = nullptr;
    } else {
      if(ce == cs->last)
	cs->last = cl;
      cl->next = ce->next;
      delete ce;
      ce = cl;
      cl = nullptr;
    }
    cs->Size--;
  }
}

void SListBase::clear()
{
  SLinkBase *cur = nullptr;
  while((cur = SListBase::get()))
    delete cur;
}

void SListBaseIter::insert(SLinkBase *ne)
{
  if(!cs->last){ // empty list
    cs->last = ne;
    ne->next = ne;
    ce = ne;
    cl = ne;
  } else {
    if((!(ce==cs->last && cl))){
      ne->next = ce->next;
      ce->next = ne;
    } else {
      ne->next = cs->last->next;
      cs->last->next = ne;
      cs->last = ne;
    }
  }
  cs->Size++;
}

SLinkBase * SListBaseIter::next() // returns next without advancing iterator
{
  if(ce){
    SLinkBase * ret = (!(ce==cs->last && cl)) ? ce->next : nullptr;
    return ret;
  }
  return nullptr;
}

SListBase const * SListBaseIter::list() const
{ return cs; }
