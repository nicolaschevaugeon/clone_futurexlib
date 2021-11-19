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
#ifndef H_SBlock
#define H_SBlock

/** A resizeable array */
template <class T>
class SBlock {
private:
  T* Elements;
  int Nelem;
public:
  ///
  SBlock( int n = 0) : Elements(new T[n]), Nelem(n) {}
  ///
  SBlock( const SBlock<T> &b);
  ///
  ~SBlock() { delete [] Elements; }
  ///
  void size(int);
  ///
  int size() const { return Nelem; }
  ///
  void reserve(int rSize);
  ///
  T& operator[](int index) { return Elements[index]; }
  ///
  T operator[](int index) const { return Elements[index]; }
};

template<class T>
SBlock<T>::SBlock( const SBlock<T> &b)
{
  Nelem = b.size();
  Elements = new T[b.size()];
  for(int i = 0; i < b.size(); i++)
    Elements[i] = b.Elements[i];
}

template<class T> 
void SBlock<T>::size(int newSize)
{
  if(newSize != size()){
    T* newEl = new T[newSize];
    if( !newEl )
      ; // error out of memory
    int nCopy = (newSize < size()) ? newSize : size();
    for(int i = 0; i < nCopy; i++)
      newEl[i] = Elements[i];
    delete [] Elements;
    Elements = newEl;
    Nelem = newSize;
  }
}

template<class T>
void SBlock<T>::reserve(int rSize)
{
  int newSize = (size() == 0) ? 1 : size();
  while( newSize < rSize)
    newSize *= 2;
  size(newSize);
}


#endif
