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

#ifndef H_Pair
#define H_Pair

/** A pair of values, the types of which can be different */
template <class L, class R>
class Pair{
private:
	L	Left;
	R	Right;
public:
	Pair() = default;
	Pair( const L& left, const R& right) : Left(left), Right(right) {}
	L	left() const { return Left; }
	void left(const L& left) { Left = left; }
	R	right() const { return Right; }
	void right(const R& right) { Right = right; }

	L	first() const { return Left; }
	R	second() const { return Right; }

};

#endif
