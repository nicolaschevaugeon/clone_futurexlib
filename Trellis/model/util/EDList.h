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
  EDList.h

  List class for embedded doubly-linked lists - items which are
  used are derived from the EDItemBase class

  Created Wed Nov 26 11:45:33 EST 1997, JDT

  $Id: EDList.h,v 1.3 2005/02/02 00:08:57 acbauer Exp $
*/

#ifndef H_EDList
#define H_EDList

class EDItemBase {
  friend class EDList;
  friend class EDListIter;
public:
  EDItemBase() { prev=nullptr; next=nullptr; }
  virtual ~EDItemBase() = default;
private:
  EDItemBase *prev;
  EDItemBase *next;
};

class EDListIter;

/** List class for embedded doubly-linked lists - items in list must
  be derived from the EDItemBase class. The list owns everything in
  it, so when the list is deleted so is everything in the list */
class EDList {
  friend class EDItemBase;
  friend class EDListIter;
public:
  EDList();
  ~EDList();
  void append(EDItemBase *item);
  void remove(EDItemBase *item);
  EDItemBase *nth(int n) const;
  inline int size() { return list_size; }
private:
  void addIter(EDListIter *iter);
  void removeIter(EDListIter *iter);

  EDItemBase *first;
  EDItemBase *last;
  int list_size;
  EDListIter **activeIters;
  int nActive;
  int nAlloc;
};

class EDListIter {
  friend class EDList;
  friend class EDItemBase;
public:
  EDListIter();
  EDListIter(const EDList *list);
  EDListIter(const EDListIter & iter);
  ~EDListIter();

  void operator=(const EDListIter &iter);
  int operator()(EDItemBase *&);
  void reset();
  EDItemBase *next();
  EDItemBase *current() const;
  EDList * list() const;
private:
  const EDList *d_list;
  EDItemBase *d_current;
};

#endif
