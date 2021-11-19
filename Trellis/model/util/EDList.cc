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
/*
  EDList.cc

  implementation of "embedded doubly-linked" lists of EDItemBase
  derived class items

  Created Wed Nov 26 11:47:32 EST 1997

  $Id: EDList.cc,v 1.3 2005/02/02 00:08:57 acbauer Exp $ 
*/

#include "EDList.h"
#include "MessageOut.h"

const int activeIterAllocSize = 10;

EDList::EDList() {

  first=nullptr;
  last=nullptr;
  list_size=0;
  activeIters=new EDListIter*[activeIterAllocSize];
  for (int i=0; i<activeIterAllocSize; i++) 
    activeIters[i]=nullptr;
  nActive=0;
  nAlloc=activeIterAllocSize;
}

EDList::~EDList() {

  if (nActive) {
    Warning("EDList::~EDList: deleting list with " << nActive << " active iterators");
  }
  // do we want to delete any iterators now or assume someone else will?
  if (nAlloc && activeIters)
    delete activeIters;
  // now delete everything in the list
  EDItemBase *item = first;
  EDItemBase *next = first;
  while(item){
    next= item->next;
    delete item;
    item = next;
  }
}

void EDList::addIter(EDListIter *iter)
{
  int i;
  // make sure there will be room to store
  if (nAlloc == nActive) {
    // make it bigger
    EDListIter **newactiveIters=new EDListIter*[2*nAlloc];
    for (i=0; i<nAlloc; i++) {
      newactiveIters[i]=activeIters[i];
      newactiveIters[i+nAlloc]=nullptr;
    }
    delete activeIters;
    activeIters=newactiveIters;
    nAlloc*=2;
  }
  // save this iter in the list of active iterators
  for (i=0; i<nAlloc; i++) {
    if (activeIters[i] == nullptr) {
      activeIters[i] = iter;
      nActive++;
      return;
    }
  }
  InternalError("EDList::addIter: should have space to store iter but don't");

}

void EDList::removeIter(EDListIter *iter)
{

  for (int i=0; i<nAlloc; i++) {
    if (activeIters[i] == iter) {
      activeIters[i]=nullptr;
      nActive--;
      return;
    }
  }
  
  Error("EDList::removeIter: iter not found in active iter list");

}

void EDList::append(EDItemBase *item) {
  EDItemBase *oldlast;

  if (!first) {
    first=item;
    last=item;
  }
  else {
    oldlast=last;
    oldlast->next=item;
    item->prev=oldlast;
    last=item;
  }

  // if there are iterators, we need to make sure they know about
  // the newly added item, if necessary
  if (nActive) {
    int i;
    for (i=0; i<nAlloc; i++) {
      EDListIter *iter=activeIters[i];
      if (iter && !iter->d_current) {
	// we appended something to the end of the list being traversed
	// by an iterator which has expired but not yet been deleted
	// need to "unexpire".  If more is appended it will be after this
	// item and we don't need any special consideration
	iter->d_current=item;
      }
    }
  }
  list_size++;
}

void EDList::remove(EDItemBase *item) {
  int i;
  EDListIter *iter;
  EDItemBase *oldprev = item->prev;
  EDItemBase *oldnext = item->next;

  // advance iterators
  if (nActive) {
    for (i=0; i<nAlloc; i++) {
      iter=activeIters[i];
      // check now if this iterator is pointing to the thing we're deleting
      if (iter && iter->d_current == item) {
	// advance the iterator, if there is a next, if not
	// move it back to the new "last"
	if (oldnext)
	  iter->d_current = oldnext;
	else
	  iter->d_current = oldprev;
      }
    }
  }

  // now just delete the item from the list
  if (oldprev)
    oldprev->next=oldnext;
  else
    first=oldnext;

  if (oldnext)
    oldnext->prev=oldprev;
  else
    last=oldprev;

  list_size--;

  // clear item's pointers for reuse in a different list
  item->prev=nullptr;
  item->next=nullptr;
 
}

EDItemBase *EDList::nth(int n) const
{
  int i;
  EDItemBase *item;

  if (n>=list_size) {
    Error("EDList::nth: requested item " << n << " from list containing " << list_size << " items");
    return nullptr;
  }

  i=0;
  item=first;
  while (i<n) {
    i++;
    item=item->next;
  }
  return item;
}

EDListIter::EDListIter()
  : d_list(nullptr), d_current(nullptr)
{}

EDListIter::EDListIter(const EDList *the_list) {
  //int i;

  d_list=the_list;
  d_current=d_list->first;

  ((EDList*)d_list)->addIter(this);
}

EDListIter::EDListIter(const EDListIter & iter)
  : d_list(iter.d_list), d_current(iter.d_current)
{
  ((EDList*)d_list)->addIter(this);
}

EDListIter::~EDListIter() 
{
  ((EDList*)d_list)->removeIter(this);
}

void EDListIter::operator=(const EDListIter &iter)
{
  if(d_list)
    ((EDList*)d_list)->removeIter(this);
  d_list = iter.d_list;
  ((EDList*)d_list)->addIter(this);
  d_current = iter.d_current;
}

int EDListIter::operator() (EDItemBase *&item) {

  item = d_current;
  if (item) {
    d_current=d_current->next;
    return 1;
  }
  return 0;
}

void EDListIter::reset() {

  d_current=d_list->first;
}

EDItemBase * EDListIter::next()
{
  EDItemBase *item = d_current;
  if(item)
    d_current = d_current->next;
  return item;
}

EDItemBase *EDListIter::current() const
{ return d_current; }

EDList * EDListIter::list() const
{ return (EDList*)d_list; }

