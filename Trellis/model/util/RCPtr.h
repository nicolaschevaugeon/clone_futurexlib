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
// classes for reference counting, pretty much taken verbatim from
// "More Effective C++" by Meyers

template<class T>
class RCPtr {
public:
  RCPtr(T* realPtr = 0);
  RCPtr(const RCPtr<T> & rhs);
  ~RCPtr();
  
  RCPtr<T> & operator = (const RCPtr<T> & rhs);
  
  T* operator->() const;
  T& operator*() const;
  
private:
  T * pointee;
  
  void init();
};

template<class T>
void RCPtr<T>::init()
{
  if(pointee==0)
    return;
  if(pointee->isShareable() == 0)
    pointee = new T(*pointee);
  pointee->addReference();
}

template<class T>
RCPtr<T>::RCPtr(const RCPtr<T> &rhs)
  : pointee(rhs.pointee)
{ init(); }

template<class T>
RCPtr<T>::RCPtr(T *realPtr)
  : pointee(realPtr)
{ init(); }

template<class T>
RCPtr<T>::~RCPtr()
{
  if(pointee)
    pointee->removeReference();
}

template<class T>
RCPtr<T> & RCPtr<T>::operator=(const RCPtr &rhs)
{
  if (pointee != rhs.pointee){
    if(pointee)
      pointee->removeReference();
    pointee = rhs.pointee;
    init();
  }
  return *this;
}

template<class T>
T* RCPtr<T>::operator->() const
{ return pointee; }

template<class T>
T& RCPtr<T>::operator*() const
{ return *pointee; }

