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

#ifndef H_Range
#define H_Range


namespace AOMD {
/** represents a range of values of the template type */
template <class T>
class Range{
private:
  T	Low;
  T	High;
public:
  Range() = default;
  Range( const T& low, const T& high) : Low(low), High(high) {}
  T	low() const { return Low; }
  void low(const T& low) { Low = low; }
  T	high() const { return High; }
  void high(const T& high) { High = high; }

  int contains(const T& value) const;
  int contains(const Range<T> & range) const;

  int operator == (const Range<T> &range) const;
};

template<class T>
int Range<T>::contains(const T& value) const
{ return ( (value >= Low) && (value <= High) ); }

template<class T>
int Range<T>::contains(const Range<T>& range) const
{ return ( (range.low() >= Low) && (range.high() <= High) ); }

template<class T>
int Range<T>::operator == (const Range<T>& range) const
{ return ( (range.low() == Low) && (range.high() == High) ); }


}
#endif
